# TSRT78 Digital signal processing

This repository is used for the labs of the course TSRT78 - Signal processing. In this course we learned:
- signal modeling of deterministic and stochastic signals;
- parametric modeling with ARMA models;
- model estimation;
- optimal filtering with Wiener and Kalman filtering;
- adaptive filtering.

Here is a link to the [course book](https://books.google.ch/books/about/Signal_Processing.html?id=UDlMygAACAAJ&source=kp_book_description&redir_esc=y), and to the [course page](https://liu.se/studieinfo/en/kurs/tsrt78/ht-2018).

## Lab 1. Speech encoding
The purpose of the lab was to investigate signal modeling in practice by modeling different signals. This was divided into three tasks were we modeled:
- a whistle
- the vowels _a_ and _o_
- a complete sentence.  
The major findings of the lab was that an AR(2)-model is good for modeling periodic signals, and that when modeling and simulating a signal, better results are achieved if the signal is divided into segments that each are modelled on their own and then composed into one simulation of the signal. The latter method is the one used in the speech encoding of Global System for Mobile Communication (GSM).

Here is a link to the [report](lab1/TSRT78_lab1_cassu286_matvi959_final.pdf).

## Lab 2. Noise canceling
This lab was about designing an adaptive filter for noise canceling. In the lab, music is corrupted by a noise signal (white noise, sum of sine-waves and chirp). 

In the first part of the lab, we model the noise signals with AR-models of different orders based on Akaikes Information Criterion. In the second part, we perform adaptive noise cancellation using LMS adaptive filter.

## Kalmar and Wiener filtering
