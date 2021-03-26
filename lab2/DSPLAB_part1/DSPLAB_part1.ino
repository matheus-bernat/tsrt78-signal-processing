/*

          DIGITAL SIGNAL PROCESSING LAB 1, PART 1
          Linköping University
          isaac.skog@liu.se
          2018-11-14
*/
#include <avr/interrupt.h>


// Buffer settings
#define DECIMATION 6      // Sample rate = 48 kHz/DECIMATION = 8kHz
#define BUFSIZE  16000    // 2 seconds of data

// AD/DA Setting
#define CANREAD() (((ADC->ADC_ISR & ADC_ISR_EOC7) && (ADC->ADC_ISR & ADC_ISR_EOC6)) ? true : false )
#define READ(chan) ( ((*(ADC->ADC_CDR + 7 - (chan))) & 0xFFF) - 0x800 )
#define READ_POT(num) ((*(ADC->ADC_CDR + 5 - num) & 0xFFF ))
#define WRITE(a, b) do { DACC->DACC_CDR = (( ((a) + 0x800) & 0xFFF) | (0x00 << 12) | (( ((b) + 0x800) & 0xFFF) << 16) | (0x01 << 28) ) ; } while ( DACC->DACC_ISR & DACC_ISR_TXRDY )


volatile int16_t y[BUFSIZE];      // y-signal (wanted signal + noise)
volatile int16_t u[BUFSIZE];      // u-signal buffer (noise source)
volatile int in_A0 = 0;       // AD channel 0
volatile int in_A1 = 0;       // AD channel 1
volatile int sample_ctr = 0;  // Sample counter used in the down sampling
volatile int buffer_ctr = 0;  // Counter used to access current buffer position
volatile boolean send_data_flag=false;

int led1 = 52;
int led2 = 50;

// ****************** SETUP  ******************* //
void setup(){

  Serial.begin(115200); // Activate the serial port

  int32_t mask_PWM_pin = digitalPinToBitMask(6);
  REG_PMC_PCER1 = 1 << 4;             // activate clock for PWM controller
  REG_PIOC_PDR |= mask_PWM_pin;       // activate peripheral functions for pin (disables all PIO functionality)
  REG_PIOC_ABSR |= mask_PWM_pin;      // choose peripheral option B
  REG_PWM_CLK = 0;                     // choose clock rate, 0 -> full MCLK as reference 84MHz
  REG_PWM_CMR7 = 0 << 9;              // select clock and polarity for PWM channel (pin7) -> (CPOL = 0)
  REG_PWM_CPRD7 = 84;                 // initialize PWM period -> T = value/84MHz (value: up to 16bit), value=8 -> 10.5MHz
  REG_PWM_CDTY7 = 42;                 // initialize duty cycle, REG_PWM_CPRD6 / value = duty cycle, for 8/4 = 50%
  REG_PWM_ENA = 1 << 7;               // enable PWM on PWM channel (pin 7 = PWML6)

  configure();  // Configure the interupts and AD/DA converters

  // Configure leds
  pinMode(led1, OUTPUT);
  digitalWrite(led1, LOW);
  pinMode(led2, OUTPUT);
  digitalWrite(led2, LOW);


  // Set pin high to power up microphone
  pinMode(2, OUTPUT);
  digitalWrite(2, HIGH);

}
// ******************************************** //


// ************** MAIN LOOP ******************* //
void loop(){
    int8_t indata=0;
    // Check if recording should start
    if (Serial.available() > 0){
      indata=Serial.read();
      if(indata=255){
        digitalWrite(led1, HIGH);
        NVIC_EnableIRQ (ADC_IRQn);
        delay(200);
        buffer_ctr = 0;
        }
      }

    // Check if data should be sent
    if(send_data_flag){

      digitalWrite(led1, LOW);
      digitalWrite(led2, HIGH);

      // Send all the data
      for (int j = 0; j < BUFSIZE; j++ ) {
      Serial.write((uint8_t)(u[j] >> 8));
      Serial.write((uint8_t)u[j]);
      Serial.write((uint8_t)(y[j] >> 8));
      Serial.write((uint8_t)y[j]);
      }

      digitalWrite(led2, LOW);
      send_data_flag=false;
     }


}
//***********************************//




//**************  SUBFUNCTIONS  ******************//

#ifdef __cplusplus
extern "C"
{
#endif
}
// AD interupt handler
void ADC_Handler (void) {

  if (CANREAD()) {
    // Add new AD reading to input variables
    in_A0 = in_A0 + READ(0);
    in_A1 = in_A1 + READ(1);


    // Increase the sample counter
    sample_ctr++;

    // Check if D samples has been collected
    if (sample_ctr == DECIMATION) {

      // Write decimated value to processing buffers
      u[buffer_ctr] = (in_A0 / DECIMATION)<<4;
      y[buffer_ctr] = (in_A1 / DECIMATION)<<4;

      // Increase the buffer counter
      buffer_ctr = buffer_ctr + 1;

      // Reset input variables and the sample counter
      in_A0 = 0;
      in_A1 = 0;
      sample_ctr = 0;

      if(buffer_ctr==BUFSIZE)
      {
      buffer_ctr=0;
      NVIC_DisableIRQ(ADC_IRQn);
      send_data_flag=true;
      }
    }
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
}

void adc_setup() {
  //NVIC_EnableIRQ (ADC_IRQn) ; // enable ADC interrupt vector
  ADC->ADC_IDR = 0xFFFFFFFF ; // disable interrupts
  ADC->ADC_IER = 0x80; // enable End-Of-Conv interrupt for AD7
  ADC->ADC_CHDR = 0xFFFF ; // disable all channels
  ADC->ADC_CHER = 0xF0 ; // enable channels AD7 and AD6
  //ADC->ADC_CGR = 0x15555555 ; // set channel gains to 1
  //ADC->ADC_CGR = 0x1555A555 ; // set channel gains to 2
  ADC->ADC_CGR = 0x1555F555 ; // set channel gains to 4
  ADC->ADC_COR = 0x000000C0 ; // set channel offsets

  ADC->ADC_MR = (ADC->ADC_MR & 0xFFFFFFF0) | (1 << 1) | ADC_MR_TRGEN | ADC_MR_ANACH ;  // 1 = trig source TIO from TC0
}
