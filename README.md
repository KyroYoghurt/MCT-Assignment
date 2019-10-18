# MCT-Assignment

## Authors

* **Naman Agarwal** - *2017A3PS0358G* - [KyroYoghurt](https://github.com/KyroYoghurt)   
* **Manish Dash** - *2017AAPS0346G* - [manishdash123](https://github.com/manishdash123)  
* **Anhadveer Khokar** - *2017A8PS0714G* - [Anhadveer](https://github.com/Anhadveer)  
* **Shirin Kaushik** - *2017AAPS0229G* - [Shirin04](https://github.com/Shirin04)  
* **Harshwardhan Shirodkar** - *2017AAPS0169G* - [harshwardhands](https://github.com/harshwardhands)  


The aim of the project is to create a GUI for data transmission over various modulation and encoding schemes and studying the BERs plots from theoretical and simulation results. 

The project is majorly divided into the following parts:

  Part: A  
    Source Data: (User Defined)  
    1. Generate Random Bits  
    2. Convert a stored text file into bits  
    3. Convert a stored gray scale image into bits 4. Convert a stored audio signal into bits  

  Part: B  
    Channel Coding (User defined)  
    1. (7,4) hamming code  
    2. BCH (127,64) code  
    3. Convolutional Code with CG=[5 7] and K=[3]  

  Part: C  
    Modulation (User defined)   
    1. 16 FSK  
    2. 16 QAM  
    3. 16 PSK  
    
  Part. D  
    Display Transmitted data parameters  
    1. Plot constellation (for 16 PSK and 16 QAM)   
    2. Plot time domain graph  
    3. Plot PSD  
    
  Part. E  
    Add noise  
    1. Add AWGN (for different SNR values)  
    
  Part F  
    Demodulate the data  
    1. Demodulate the symbols and convert into bits  
    
  Part G  
    Display the results  
    1. Plot BER curves and Compare the theoretical and Simulated Results with and without coding.  
    2. Reconstruct the Text, Image and Audio for different SNR values and display/ play it to the user (with and without              coding)  

## Getting Started

Download the .zip file and extract the content. The unzipped folder will contain:  
1. Main .m file: gui_3.m  
2. GUI figure: fig_3.fig  
3. Sample image file: 
4. Sample text file:
5. Sample audio file    
To test the code, downlaod the attached text, image and audio inputs as well. Ensure that all the files are at a common location and the path to the location has been added to the working folder of MATLAB. 

### Prerequisites

Latest version of MATLAB with the following toolboxes:  
1. Communication toolbox  
2. Digital Signal Processing toolbox  
3. Image Processing toolbox  


## Running the tests

1. Open the .m file and RUN it.  
2. A GUI will pop up.  
3. Select one of the 4 options in 'Source Data': Generate Random Bits, Text file to bits, Grayscale image to bits, Audio to bits. The input is imported from the downloaded sample inputs and converted to binary data.  
4. Next, select a 'Channel Coding' scheme: (7,4) Hamming Code, BCH(127,64) code or Convolutional Code. Do not select multiple options or same option multiple times. The processing may take some time (3-5 minutes) depending on the length of input signal and processor.  
5. Select a modulation scheme: 16 FSK, 16 PSK, 16 QAM and wait for sometime. The BER plot for the selected modulation scheme will appear comparing the theoretical and simulation results for uncoded and coded data. 16 QAM will plot the fastest. Two more plots displaying the time domain graph of the input signal and Power Spectral Density(PSD) will appear. For 16 PSK and 16 QAM, a constellation plot can be observed.  
6. In the decoding section, select a decoding scheme. Ensure that the decoding scheme is same as the previously chosen coding scheme. For example, if the encoding scheme is (7,4) Hamming Code, the decoding scheme should also be (7,4) Hamming Code.  
7. Once the decoding is complete, reconstruct the decoded signal by selecting the input signal type in the 'Reconstruction' section. The reconstructed input will be displayed as following if the input was:  
    a. 'Random bits', an array of bits will appear in command window  
    b. 'Text file', a string will appear in command window  
    c. 'Grayscale Image', an image will appear in the GUI  
    d. 'Audio', an audio output will be heard from the system's inbuilt speakers.  


## Deployment

The project can be deployed into a live system by altering the input sources. The sample text file can be upated dynamically for continuous text input. A new file path can be added in the callback function for image in 'Source Data' section. Similarly, audioread() function can take input from the device's microphone and transmit through the selected channel. 

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgment

* Mr. Naveen Gupta, Asst. Professor (Grade-I), Instructor for 'Modern Communication Technologies'(ECE F418),  BITS Pilani, KK Birla Goa Campus
