/*
 *
 * Fixed point implmentation of a leaky LMS algorithm for noise suppression.
 *
 * Isaac Skog, isaac.skog@liu.se, 2017-12-25
 *
 */
#include <avr/interrupt.h>

// Filter settings
#define NR_OF_FILTER_TAPS 40      // Equivalent to n_b in the lab instructions
#define DELAY 9                   // Equivalent to n_k in the lab instructions
//#define gamma 32767                   // Equivalent to n_k in the lab instructions


// Buffer settings
#define DECIMATION 6      // Sample rate = 48 kHz/DECIMATION
#define BUFSIZE 128       // MAX_NR_TAPS+MAX_DELAY
#define BUFMASK 127       // Mask used in the circular buffer
#define MAX_NR_TAPS 64    // Max nb
#define MAX_DELAY 64      // Max nk
#define MAX_INT16 32767
#define MIN_INT16 -32768


#if (NR_OF_FILTER_TAPS >MAX_NR_TAPS) | (NR_OF_FILTER_TAPS < 1)
#error Not a valid number of filter taps! Most be between 1 and 64
#endif

#if (DELAY > MAX_DELAY) | (DELAY < 0)
#error Not a valid delay! The delay most be between 0 and 64
#endif

// AD/DA Setting
#define CANREAD() (((ADC->ADC_ISR & ADC_ISR_EOC7) && (ADC->ADC_ISR & ADC_ISR_EOC6)) ? true : false )
#define READ(chan) ( ((*(ADC->ADC_CDR + 7 - (chan))) & 0xFFF) - 0x800 )
#define READ_POT(num) ((*(ADC->ADC_CDR + 5 - num) & 0xFFF ))
#define WRITE(a, b) do { DACC->DACC_CDR = (( ((a) + 0x800) & 0xFFF) | (0x00 << 12) | (( ((b) + 0x800) & 0xFFF) << 16) | (0x01 << 28) ) ; } while ( DACC->DACC_ISR & DACC_ISR_TXRDY )


volatile int u[2 * BUFSIZE];        // u-signal buffer (noise source)
volatile int y;                     // y-signal (utility signal + noise)
volatile int s_hat;                 // Estimated noise free signal
volatile int in_A0 = 0;             // AD channel 0
volatile int in_A1 = 0;             // AD channel 1
volatile int out_A0 = 0;            // DA channel 0
volatile int out_A1 = 0;            // DA channel 1
volatile int pot0 = 0;              // Potentiometer value of left potentiometer
volatile int pot1 = 0;              // Potentiometer value of right potentiometer
volatile int sample_ctr = 0;        // Sample counter used in the down sampling
volatile int buffer_ctr = 0;        // Counter used to access current buffer position
volatile int w[MAX_NR_TAPS];        // Vector with the filter parameters
volatile uint8_t mu_shift = 6;      // Step size 1/(2^sb)
volatile int lambda=MAX_INT16;      // Leakage factor
int buttonPins [] = {34, 36} ;

boolean lastButtons [] = {true, true} ;
boolean currentButtons [] = {true, true} ;
boolean process_data_flag = false;
volatile int state = 0;
int lastTime = 0;
int ledPins [] = {52, 50, 48, 46, 44, 42, 40, 38} ;



// ****************** SETUP  ******************* //
void setup(){

  int32_t mask_PWM_pin = digitalPinToBitMask(6);
  REG_PMC_PCER1 = 1 << 4;             // activate clock for PWM controller
  REG_PIOC_PDR |= mask_PWM_pin;       // activate peripheral functions for pin (disables all PIO functionality)
  REG_PIOC_ABSR |= mask_PWM_pin;      // choose peripheral option B
  REG_PWM_CLK = 0;                     // choose clock rate, 0 -> full MCLK as reference 84MHz
  REG_PWM_CMR7 = 0 << 9;              // select clock and polarity for PWM channel (pin7) -> (CPOL = 0)
  REG_PWM_CPRD7 = 84;                 // initialize PWM period -> T = value/84MHz (value: up to 16bit), value=8 -> 10.5MHz
  REG_PWM_CDTY7 = 42;                 // initialize duty cycle, REG_PWM_CPRD6 / value = duty cycle, for 8/4 = 50%
  REG_PWM_ENA = 1 << 7;               // enable PWM on PWM channel (pin 7 = PWML6)


  for (int j = 0; j < 8; j++ ) {      // Make sure leds are turned off
    pinMode (ledPins [j], OUTPUT) ;
    digitalWrite(ledPins[j], LOW);
  }

  // Flush buffers and filter parameters
  for (int j = 0; j < 2 * BUFSIZE; j++) {
    u[j] = 0;
  }

  // Flush filter parameters
  for (int j = 0; j < MAX_NR_TAPS; j++) {
    w[j] = 0;
  }

   for (int j = 0; j < 2; j++ ) {
    pinMode (buttonPins [j], INPUT) ;
    digitalWrite (buttonPins [j], HIGH) ;
    attachInterrupt (digitalPinToInterrupt(buttonPins [j]), buttonPressed, CHANGE) ;
  }


  // Configure the interupts and AD/DA converters
  configure();

  // Display the current operation state/mode
  digitalWrite(ledPins[state], HIGH);

  // Set pin high to power up microphone
  pinMode(2, OUTPUT);
  digitalWrite(2, HIGH);
}
// ******************************************** //


// ************** MAIN LOOP ******************* //
void loop(){


  /* LMS algorithm
   *
   *  The filter taps are assumed to be in the range [-1,1] and are represented using 16 bits.
   *
  */
  if ((process_data_flag == HIGH) ) {
    process_data_flag = LOW;
    //digitalWrite(2, HIGH);

    // Estimate noise
    int y_hat = 0;
    for (int j = 0; j < NR_OF_FILTER_TAPS; j++ ) {
      y_hat = y_hat  + ((w[j] * u[buffer_ctr + BUFSIZE - j - DELAY])  >> 15);
    }

    // Estimate the noise free signal
    s_hat = (y - y_hat);

    // Update the filter parameters. Note that the step length is implemented as a down shift.
    for (int j = 0; j < NR_OF_FILTER_TAPS; j++ ){

      w[j] = (lambda*w[j] + ((u[buffer_ctr + BUFSIZE - j - DELAY] * s_hat) >> mu_shift) ) >> 15;

      // Check that the parameters are within the valid range, else truncate to maximum/minimum value
      if (w[j] > MAX_INT16)
      {
        w[j] = MAX_INT16;
      }
      if (w[j] < MIN_INT16)
      {
        w[j] = MIN_INT16;
      }
    }

    // Check which signal that should be outputted
    if (state == 0) {
      // Orignal sound
      out_A0 = y >> 4;  // Down shift to 12 bits (DA resolution)
      out_A1 = y >> 4;
    }
    else if (state == 1)
    {
      // Estimated "noise" free signal
      out_A0 = s_hat >> 4;
      out_A1 = s_hat >> 4;
    }
    else
    {
      // Estimated noise
      out_A0 = y_hat >> 4;
      out_A1 = y_hat >> 4;
    }
    //digitalWrite(2, LOW);
  }


  // Read potentiometer to adjust the step length and leakage factor
  pot0 = READ_POT(0);
  pot1 = READ_POT(1);

  // Find the corresponding step length
  int level = 820;
  // Turn off led
  digitalWrite(ledPins[mu_shift], LOW);
  for (uint8_t j = 3; j < 8; j++ ) {
    if (pot0 < level) {
      mu_shift = j;
      // Turn on led that displace the current step length
      digitalWrite(ledPins[mu_shift], HIGH);
      break;
    }
    level = level+820;
  }

  // Leakage factor ( 0- 2^10/2^12)
  lambda=MAX_INT16-(pot1>>2);

}
//***********************************//




//**************  SUBFUNCTIONS  ******************//

#ifdef __cplusplus
extern "C"
{
#endif
}

void buttonPressed(){
  currentButtons[0] = digitalRead (buttonPins[0]) ;
  currentButtons[1] = digitalRead (buttonPins[1]) ;

  if ( ((lastButtons[0] & !currentButtons[0]) ^ (lastButtons[1] & !currentButtons[1])) && ((millis() - lastTime) > 50)) {

    // Turn of led
    digitalWrite(ledPins[state], LOW);

    if (currentButtons[1]) {
      if (state != 2) state ++ ;
      else state = 0;
    } else if (currentButtons[0]) {
      if (state != 0) state-- ;
      else state = 2;
    }

    // Turn on led
    digitalWrite(ledPins[state], HIGH);
  }
  lastTime = millis();
  lastButtons[0] = currentButtons[0] ;
  lastButtons[1] = currentButtons[1] ;
}


// AD interupt handler
void ADC_Handler (void) {

  if (CANREAD()) {
    //noInterrupts();
    // Add new AD reading to input variables
    in_A0 = in_A0 + READ(0);
    in_A1 = in_A1 + READ(1);


    // Increase the sample counter
    sample_ctr++;

    // Check if D samples has been collected
    if (sample_ctr == DECIMATION) {

      // Write output to DAC
      WRITE(out_A0, out_A1);

      // Increase the buffer counter
      buffer_ctr = (buffer_ctr + 1) & BUFMASK;

      // Write decimated value to processing buffers
      u[buffer_ctr] = (in_A0 / DECIMATION) << 4;
      u[buffer_ctr + BUFSIZE] = u[buffer_ctr];
      y = (in_A1 / DECIMATION) << 4;

      // Reset input variables and the sample counter
      in_A0 = 0;
      in_A1 = 0;
      sample_ctr = 0;

      // Signal that data should be processed
      process_data_flag = HIGH;
    }
    //interrupts();
  }
}

void configure() {
  adc_setup() ;
  pmc_enable_periph_clk (TC_INTERFACE_ID + 0 * 3 + 0) ; //clock the TC0 channel 0
  TcChannel * t = &(TC0->TC_CHANNEL)[0] ;  // pointer to TC0 registers for its channel 0
  t->TC_CCR = TC_CCR_CLKDIS ;  // disable internal clocking while setup regs
  t->TC_IDR = 0xFFFFFFFF ;  // disable interrupts
  t->TC_SR ;  // read int status reg to clear pending
  t->TC_CMR = TC_CMR_TCCLKS_TIMER_CLOCK1 |  // use TCLK1 (prescale by 2, = 42MHz)
              TC_CMR_WAVE |  // waveform mode
              TC_CMR_WAVSEL_UP_RC |  // count-up PWM using RC as threshold
              TC_CMR_EEVT_XC0 |  // Set external events from XC0 (this setup TIOB as output)
              TC_CMR_ACPA_CLEAR | TC_CMR_ACPC_CLEAR |
              TC_CMR_BCPB_CLEAR | TC_CMR_BCPC_CLEAR ;

  t->TC_RC =  875 ;  // counter resets on RC, so sets period in terms of 42MHz clock (48kHz sampling frequency)
  t->TC_RA =  440 ;  // roughly square wave
  t->TC_CMR = (t->TC_CMR & 0xFFF0FFFF) | TC_CMR_ACPA_CLEAR | TC_CMR_ACPC_SET ;  // set clear and set from RA and RC compares
  t->TC_CCR = TC_CCR_CLKEN | TC_CCR_SWTRG ;  // re-enable local clocking and switch to hardware trigger source
  dac_setup() ;
}

void adc_setup() {
  NVIC_EnableIRQ (ADC_IRQn) ; // enable ADC interrupt vector
  ADC->ADC_IDR = 0xFFFFFFFF ; // disable interrupts
  ADC->ADC_IER = 0x80 ; // enable End-Of-Conv interrupt for AD7
  ADC->ADC_CHDR = 0xFFFF ; // disable all channels
  ADC->ADC_CHER = 0xF0 ; // enable channels AD7 and AD6
  //ADC->ADC_CGR = 0x15555555 ; // set channel gains to 1
  //ADC->ADC_CGR = 0x1555A555 ; // set channel gains to 2
  ADC->ADC_CGR = 0x1555F555 ; // set channel gains to 4
  ADC->ADC_COR = 0x000000C0 ; // set channel offsets

  ADC->ADC_MR = (ADC->ADC_MR & 0xFFFFFFF0) | (1 << 1) | ADC_MR_TRGEN | ADC_MR_ANACH ;  // 1 = trig source TIO from TC0

}

void dac_setup() {
  pmc_enable_periph_clk (DACC_INTERFACE_ID) ; // start clocking DAC
  DACC->DACC_CR = DACC_CR_SWRST ; // reset DAC
  DACC->DACC_MR = DACC_MR_TAG | DACC_MR_WORD | DACC_MR_REFRESH (0x0F) | (24 << DACC_MR_STARTUP_Pos) ; // configure DACC Mode Register
  DACC->DACC_IDR = 0xFFFFFFFF ; // disable all interrupts
  DACC->DACC_CHER = 0x3 ; // enable DAC0 and DAC1
}
