# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 20:28:06 2020

@author: Arthur Cruz Morbach
"""

import matplotlib.pyplot as plt
import numpy as np
from math import pi
import wave as fw
import wavio

def closew():
    plt.close("all")

closew()

def fftfun(x,Fs):
    n = len(x) #tamanho do sinal
    k = np.arange(n) #vetor em k
    T = n/Fs
    frq = k/T #os dois lados do vetor de frequência
    frq = frq[range(int(n/2))] #apenas um lado
    
    X = 2*np.fft.fft(x)/n #cálculo da fft e normalização por n
    xl = X[range(int(n/2))]
    return(frq, abs(xl))


#----------------------------


def criaWAV(out_filename, data, Damp):
    #out_file = xxxxx.wav    
    data = data*Damp
    wavio.write(out_filename, data, 44100, sampwidth = 2)#sampwidth=arquivoWav.getsampwidth())
    return('ok file')


fs = 44100  #taxa de amostragem


#---------------Low Frequency Shelving Filter--------------------------------------------------
def LowFreqShelvingBoost(x,fc,G):
    #x = sinal de entrada    
    #fc = frequência de corte
    #G = Ganho em dB
    K = np.tan((pi*fc)/fs)
    V0 = 10**(G/20)
    b0 = (1+((2*V0)**0.5)*K+V0*(K**2))/(1+(2**0.5)*K+K**2)
    b1 = (2*((V0*K**2)-1))/(1+(2**0.5)*K+K**2)
    b2 = (1-((2*V0)**0.5)*K+V0*K**2)/(1+(2**0.5)*K+K**2)
    a1 = (2*((K**2)-1))/(1+(2**0.5)*K+K**2)
    a2 = (1-(2**0.5)*K+K**2)/(1+(2**0.5)*K+K**2)
    
    #Variáveis de estado
    xh1=0;
    xh2=0;
        
    yh1=0;
    yh2=0;
        
    y = np.zeros(len(x));
        
    #Convolução Filtro
    for n in range(0,len(x)):
        y[n]= b0*x[n] + b1*xh1 + b2*xh2 - a1*yh1 - a2*yh2;
            
        yh2=yh1;
        yh1=y[n];
            
        xh2=xh1;
        xh1=x[n];
    return(y)
    
def LowFreqShelvingCut(x,fc,G):
    #x = sinal de entrada    
    #fc = frequência de corte
    #G = Ganho em dB
    K = np.tan((pi*fc)/fs)
    V0 = 10**(G/20)
    b0 = (1+(2**0.5)*K+K**2)/(1+((2*V0)**0.5)*K+(V0*K**2))
    b1 = (2*((K**2)-1))/(1+((2*V0)**0.5)*K+(V0*K**2))
    b2 = (1-(2**0.5)*K+K**2)/(1+((2*V0)**0.5)*K+(V0*K**2))
    a1 = (2*((V0*K**2)-1))/(1+((2*V0)**0.5)*K+(V0*K**2))
    a2 = (1-((2*V0)**0.5)*K+V0*K**2)/(1+((2*V0)**0.5)*K+(V0*K**2))
    
    #Variáveis de estado
    xh1=0;
    xh2=0;
        
    yh1=0;
    yh2=0;
        
    y = np.zeros(len(x));
        
    #Convolução Filtro
    for n in range(0,len(x)):
        y[n]= b0*x[n] + b1*xh1 + b2*xh2 - a1*yh1 - a2*yh2;
            
        yh2=yh1;
        yh1=y[n];
            
        xh2=xh1;
        xh1=x[n];
    return(y)
    
    
#-------------BAND PASS FILTER--------------------------------------------
def PeakBoost(x,fc,G):
    #x = sinal de entrada    
    #fc = frequência de corte
    #G = Ganho em dB
    Q = 15
    K = np.tan((pi*fc)/fs)
    V0 = 10**(G/20)
    b0 = (1+(V0/Q)*K+K**2)/(1+(1/Q)*K+K**2)
    b1 = (2*(K**2-1))/(1+(1/Q)*K+K**2)
    b2 = (1-(V0/Q)*K+K**2)/(1+(1/Q)*K+K**2)
    a1 = (2*(K**2-1))/(1+(1/Q)*K+K**2)
    a2 = (1-(1/Q)*K+K**2)/(1+(1/Q)*K+K**2)
    
    #Variáveis de estado
    xh1=0;
    xh2=0;
        
    yh1=0;
    yh2=0;
        
    y = np.zeros(len(x));
        
    #Convolução Filtro
    for n in range(0,len(x)):
        y[n]= b0*x[n] + b1*xh1 + b2*xh2 - a1*yh1 - a2*yh2;
            
        yh2=yh1;
        yh1=y[n];
            
        xh2=xh1;
        xh1=x[n];
    return(y)
    
def PeakCut(x,fc,G):
    #x = sinal de entrada    
    #fc = frequência de corte
    #G = Ganho em dB
    Q = 15
    K = np.tan((pi*fc)/fs)
    V0 = 10**(G/20)
    b0 = (1+(1/Q)*K+K**2)/(1+(V0/Q)*K+K**2)
    b1 = (2*(K**2-1))/(1+(V0/Q)*K+K**2)
    b2 = (1-(1/Q)*K+K**2)/(1+(V0/Q)*K+K**2)
    a1 = (2*(K**2-1))/(1+(V0/Q)*K+K**2)
    a2 = (1-(V0/Q)*K+K**2)/(1+(V0/Q)*K+K**2)
    
    #Variáveis de estado
    xh1=0;
    xh2=0;
        
    yh1=0;
    yh2=0;
        
    y = np.zeros(len(x));
        
    #Convolução Filtro
    for n in range(0,len(x)):
        y[n]= b0*x[n] + b1*xh1 + b2*xh2 - a1*yh1 - a2*yh2;
            
        yh2=yh1;
        yh1=y[n];
            
        xh2=xh1;
        xh1=x[n];
    return(y)
    
    
#---------------High Frequency Shelving Filter--------------------------------------------------
def HighFreqShelvingBoost(x,fc,G):
    #x = sinal de entrada    
    #fc = frequência de corte
    #G = Ganho em dB
    K = np.tan((pi*fc)/fs)
    V0 = 10**(G/20)
    b0 = (V0+((2*V0)**0.5)*K+K**2)/(1+(2**0.5)*K+K**2)
    b1 = (2*((K**2)-V0))/(1+(2**0.5)*K+K**2)
    b2 = (V0-((2*V0)**0.5)*K+K**2)/(1+(2**0.5)*K+K**2)
    a1 = (2*((K**2)-1))/(1+(2**0.5)*K+K**2)
    a2 = (1-(2**0.5)*K+K**2)/(1+(2**0.5)*K+K**2)
    
    #Variáveis de estado
    xh1=0;
    xh2=0;
        
    yh1=0;
    yh2=0;
        
    y = np.zeros(len(x));
        
#Convolução Filtro
    for n in range(0,len(x)):
        y[n]= b0*x[n] + b1*xh1 + b2*xh2 - a1*yh1 - a2*yh2;
            
        yh2=yh1;
        yh1=y[n];
            
        xh2=xh1;
        xh1=x[n];
    return(y)
    
    
def HighFreqShelvingCut(x,fc,G):
    #x = sinal de entrada    
    #fc = frequência de corte
    #G = Ganho em dB
    K = np.tan((pi*fc)/fs)
    V0 = 10**(G/20)
    b0 = (1+(2**0.5)*K+K**2)/(V0+((2*V0)**0.5)*K+K**2)
    b1 = (2*((K**2)-1))/(V0+((2*V0)**0.5)*K+K**2)
    b2 = (1-(2**0.5)*K+K**2)/(V0+((2*V0)**0.5)*K+K**2)
    a1 = (2*(((K**2)/V0)-1))/(1+((2/V0)**0.5)*K+((K**2)/V0))
    a2 = (1-((2/V0)**0.5)*K+((K**2)/V0))/(1+((2/V0)**0.5)*K+((K**2)/V0))
    
    #Variáveis de estado
    xh1=0;
    xh2=0;
        
    yh1=0;
    yh2=0;
        
    y = np.zeros(len(x));
        
#Convolução Filtro
    for n in range(0,len(x)):
        y[n]= b0*x[n] + b1*xh1 + b2*xh2 - a1*yh1 - a2*yh2;
            
        yh2=yh1;
        yh1=y[n];
            
        xh2=xh1;
        xh1=x[n];
    return(y)


def equalizador(filename, out_filename, LowFilter_G, PeakFilter_G, PeakFilter_G1, PeakFilter_G2, PeakFilter_G3, PeakFilter_G4, PeakFilter_G5, PeakFilter_G6, PeakFilter_G7, PeakFilter_G8, PeakFilter_G9, PeakFilter_G10, PeakFilter_G11,HighFilter_G):
    #filename = 'C:/Users/xxxx/xxxx/arquivo_de_entrada.wav'
    #out_filename = 'file_name.wav'
    #------------Abrindo arquivo de áudio-----------------------------------------

    arqff = filename;
    arquivoWav = fw.open(arqff, 'r');
    
    print("Número canais: ", arquivoWav.getnchannels());
    print("Número bytes: ", arquivoWav.getsampwidth());
    print("Taxa de amostragem: ", arquivoWav.getframerate());
    print("Número de frames: ", arquivoWav.getnframes());
    print("Compactação: ", arquivoWav.getcompname());
    
    tipos = np.uint8;
    Damp = 256;
    if arquivoWav.getsampwidth()>1:
        tipos = np.int16;
        Damp = 32760;
    
    frames = arquivoWav.readframes(-1);
    x = np.fromstring(frames, tipos)/Damp;
    fim = np.size(x);
    
    fs = 44100  #arquivo.Wav.getframerate()  #taxa de amostragem
    
    dt = 1.0 / arquivoWav.getframerate();
    t = np.arange(0, arquivoWav.getnframes() * dt, dt);
    arquivoWav.close();
    
    
#-------Processo-------------------------------------------------
    
    if(LowFilter_G>=0):
        y_l = LowFreqShelvingBoost(x, 100, abs(LowFilter_G))
    else:
        y_l = LowFreqShelvingCut(x, 100, abs(LowFilter_G))
        
        
    if(PeakFilter_G>=0):
        y_p = LowFreqShelvingBoost(y_l, 400, abs(PeakFilter_G))
    else:
        y_p = LowFreqShelvingCut(y_l, 400, abs(PeakFilter_G))    
        
        
    if(PeakFilter_G1>=0):
        y_p1 = PeakBoost(y_p, 600, abs(PeakFilter_G1))
    else:
        y_p1 = PeakCut(y_p, 600, abs(PeakFilter_G1))
        
        
    if(PeakFilter_G2>=0):
        y_p2 = PeakBoost(y_p1, 1e3, abs(PeakFilter_G2))
    else:
        y_p2 = PeakCut(y_p1, 1e3, abs(PeakFilter_G2))
        
        
    if(PeakFilter_G3>=0):
        y_p3 = PeakBoost(y_p2, 2e3, abs(PeakFilter_G3))
    else:
        y_p3 = PeakCut(y_p2, 2e3, abs(PeakFilter_G3))
        
        
    if(PeakFilter_G4>=0):
        y_p4 = PeakBoost(y_p3, 3e3, abs(PeakFilter_G4))
    else:
        y_p4 = PeakCut(y_p3, 3e3, abs(PeakFilter_G4))
        
        
    if(PeakFilter_G5>=0):
        y_p5 = PeakBoost(y_p4, 4e3, abs(PeakFilter_G5))
    else:
        y_p5 = PeakCut(y_p4, 4e3, abs(PeakFilter_G5))
        
        
    if(PeakFilter_G6>=0):
        y_p6 = PeakBoost(y_p5, 5e3, abs(PeakFilter_G6))
    else:
        y_p6 = PeakCut(y_p5, 5e3, abs(PeakFilter_G6))
        
        
    if(PeakFilter_G7>=0):
        y_p7 = PeakBoost(y_p6, 8e3, abs(PeakFilter_G7))
    else:
        y_p7 = PeakCut(y_p6, 8e3, abs(PeakFilter_G7))
        
        
    if(PeakFilter_G8>=0):
        y_p8 = PeakBoost(y_p7, 10e3, abs(PeakFilter_G8))
    else:
        y_p8 = PeakCut(y_p7, 10e3, abs(PeakFilter_G8))
        
        
    if(PeakFilter_G9>=0):
        y_p9 = PeakBoost(y_p8, 12e3, abs(PeakFilter_G9))
    else:
        y_p9 = PeakCut(y_p8, 12e3, abs(PeakFilter_G9))
        
        
    if(PeakFilter_G10>=0):
        y_p10 = PeakBoost(y_p9, 15e3, abs(PeakFilter_G10))
    else:
        y_p10 = PeakCut(y_p9, 15e3, abs(PeakFilter_G10))
        
        
    if(PeakFilter_G11>=0):
        y_p11 = PeakBoost(y_p10, 17e3, abs(PeakFilter_G11))
    else:
        y_p11 = PeakCut(y_p10, 17e3, abs(PeakFilter_G11))
               
        
    if(HighFilter_G>=0):
        y = HighFreqShelvingBoost(y_p11, 19e3, abs(HighFilter_G))
    else:
        y = HighFreqShelvingCut(y_p11, 19e3, abs(HighFilter_G)) 
    
    #-----FFT
    
    frq, X = fftfun(x,fs)
    X = 20*np.log10(X)

    frq1, Y = fftfun(y,fs)
    Y = 20*np.log10(Y)
        
    #---------------PLOTS------------------------
    
    plt.figure(0) 
    plt.subplot(2,2,1)
    plt.plot(t,x)
    plt.title('Signal')
    plt.xlabel('t (s)')
    plt.ylabel('Amplitude')
    plt.grid('on')
    
    plt.subplot(2,2,3)
    plt.plot(frq,X,'r')
    plt.title('Fourier - Signal')
    plt.xlabel('Freq (Hz)')
    plt.ylabel('|X(f)|dB')
    plt.grid('on')
    
    plt.subplot(2,2,2)
    plt.plot(t,y)
    plt.title('Filtered Signal')
    plt.xlabel('t (s)')
    plt.ylabel('Amplitude')
    plt.grid('on')
    
    plt.subplot(2,2,4)
    plt.plot(frq1,Y,'r')
    plt.title('Fourier - Filtered Signal')
    plt.xlabel('Freq (Hz)')
    plt.ylabel('|X(f)|dB')
    plt.grid('on')
    
    criaWAV(out_filename, y, Damp)
    
    return(1)
    
    
