# Gabor-Transform

**Abstract**

As we demonstrate in the Ultrasound example, Fourier transform is very useful for analyzing frequency information for stationary data. However, when it comes to non-stationary data, i.e. the frequency evolves along time, Fourier transform loses all the time information. Thus, Gabor Transform (GT), the short-window Fourier Transform is created. Instead of taking Fourier Transform to the whole data at once, GT slices up the time, focuses on certain time-window data, and do Fourier Transform on this window, to gain information on frequency. By doing this way, we are able to extract time information and frequency information, i.e. we are able to tell which frequencies happen during which time window. We use three simple music pieces to demonstrate how to use GT.

**Introduction and Overview**

Here, we use three music pieces to illustrate how Gabor Transform can be used to extract both time and frequency information.

**Problem Description**

Typically, when a music piece is given to us, we only have the amplitude information for a certain period of time. We don't have information on which frequency is played at which time point, which is what we need if we hope to reconstruct the whole piece. As we mentioned above, GT can help us with this problem. However, a drawback of GT is that there is a trade off between the time information and frequency information. If we use very small time window, we are able to gain more accurate time for each frequency to be played at, but we will lose low frequency information since they typically need large time window to be detected. Another problem when analyzing music piece is that instruments, except recorder, would produce overtone which makes it more difficult to get the accurate information on frequency.

**General Approach**

**Handel Example**

A simple Garbor transform would be applied to the first piece to demonstrate how to gain the spectrogram of a sound piece. Then, the window width of the Gabor filter would be adjusted to illustrate the trade off between time and frequency information. Afterwards, the translation, i.e. the length of time steps we use (if the translation is small, we would sample more data since we move slowly from the beginning to the end point; if the translation is large, we would sample less data, since we are making big jumps between the beginning and the ending point), would be adjusted to illustrate how over-sampling and under-sampling work. Different Gabor windows are used here to demonstrate their different uses: Gaussian window, Mexican hat wavelet, step-function, Shannon window.

**Piano and Recorder**
"Mary had a little lamb" would be played on both piano and recorder. To see the difference between piano and recorder, two spectrograms are produced. We can see that compared to piano, recorder does not have the overtone structure, meaning there is no layers of frequencies appearing in the spectrogram. To obtain a good music score for the piano piece, GT will first be used to obtain the spectrogram for us to inspect which frequencies are producing overtones, then Fourier transform is used to filter out frequencies overtones.

**Theoretical Background**

As we learned from our textbook \cite{kutz_2013}, to solve the loss of time information problem, Gabor introduced the way to slice up time, use Gabor filter around certain time window, and then apply Fourier transform. Gabor transform, or short-time Fourier Transform is defined as follows:

    G[f](t,\omega) = \Tilde{f}_g(t, \omega)=\int^{\infty}_{-\infty}f(\tau)g(\tau-t)e^{-i\omega\tau}d\tau(f,\bar{g}_{t,\omega})

This would be our main tool for this whole project. A main Gabor filter we would be using here is Gaussian window:

    g(t)=e^{-a(t-b)^2}

where a is the window width and b is the time point we can change to slide through the whole time duration, which we can adjust to extact time and frequency information.

To solve the drawback of GT, wavelet is introduced. The basic idea of wavelet is that we don't have to decide which a, or window width, we want to use when filtering. When using wavelet, a certain range of the window width would be used to extract more accurate information, thus resolving the trade off issue we would have when using Garbor transform. As illustrated in our textbook \cite{kutz_2013}, one good example is Mexican Hat wavelet:

    \psi(t)=(1-t^2)e^{-t^2/2}=-d^2/dt^2(e^{-t^2/2})

    \hat{\psi}(\omega)=\hat{\psi}_{1,0}(\omega)=\sqrt{2\pi}\omega^2e^{-\omega^2/2}


**Algorithm Implementation and Development**

Here, the way the code is developed is illustrated.

handel music is loaded from Matlab to y, we then transpose it as store in v, which is a necessary step for applying Gabor filter later. The final point was taken out since we want to make sure the length is even for faster fft.

the music time length, time domain and frequency domain is created.

The following illustrate how Gabor filter is applied to the music piece to obtain spectrogram.
    

    \STATE{Define \texttt{tslide} as a time vector from 0 to the time length of the music piece with distance $\Delta t$}
    \STATE{Define $vgt spec$ to store the Fourier transformed and filtered result.}
    \FOR{$j = 1:length(tslide)$}
        \STATE{Define $g$ as the Gaussian window for filtering out information outside of certain window width}
        \STATE{Multiply $g$ with v and store the result as $vg$}
        \STATE{Fourier Transform $vg$, $vgt = fft(vg)$}
        \STATE{Add the result to $vgt spec$.}
    \ENDFOR
    \STATE{Visualize $vgt spec$, which would give us the spectrogram.}

    
The window width is then changed to $50$ and $.5$ in Algorithm~\ref{alg:part1_1} to see the difference of the spectrograms and analyze the trade off.

The distance $\Delta t$ is then changed to $1.0$ and $.05$ in Algorithm~\ref{alg:part1_1} to see the difference of under-sampling and over-sampling.

The filter function is then changed to illustrate different Gabor windows in Algorithm~\ref{alg:part1_1}, such as Gaussian window as in Algorithm~\ref{alg:part1_1}, Mexican hat wavelet, Shannon wavelet, step function. Different filter function $g$ is provided here.


    g_1 = e^{-b(t-tslide(j)^2)}

    g_2=(1-(\frac{t-tslide(j)}{b})^2)\frac{2}{\sqrt{3b}\pi^{\frac{1}{4}}}e^{-\frac{(t-tslide(j))^2}{2b^2}}

    g_3=sinc(\frac{t-tslide(j)}{2})cos(3\pi\frac{t-tslide(j)}{2})

    g_4= \begin{cases} 
      1 & tslide(j) - \Delta t \leq t\leq tslide(j) \\
      0 & otherwise \\

The following would be for analyzing the piano and recorder cases. We first use Algorithm~\ref{alg:part1_1} above to create spectrograms for both pieces and observe their differences.

Then, Gaussian filter is applied after Fourier Transform to filter out overtones using Algorithm~\ref{alg:part2_1}. We can get a very clean result showing us the music notes played during certain time period.
    

    \STATE{Define \texttt{tslide} as a time vector from 0 to the time length of the music piece with distance $\Delta t$}
    \STATE{Define $vgt spec$ to store the Fourier transformed and filtered result.}
    \STATE{Define \texttt{omegas} to store the frequency of the spike value in $vgt spec$}
    \STATE{Define a vector storing different window width we want to try to get different information.}
    \FOR{$i=1:length(width vector)$}
    \FOR{$j = 1:length(tslide)$}
        \STATE{Define $g$ as the Gaussian window for filtering out information outside of certain window width}
        \STATE{Multiply $g$ with v and store the result as $vg$}
        \STATE{Fourier Transform $vg$, $vgt = fft(vg)$}
        \STATE{Add the result to $vgt spec$.}
        \STATE{Find the frequency, $\omega$, of the spike value in $vgt$}
        \STATE{Use the $\omega$ found above to define the Gaussian filter $f$}
        \STATE{Multiply $f$ with $vgt$ and store the result as $vgft$}
        \STATE{Visualize the filtered main frequency}
    \ENDFOR
    \ENDFOR

**Computational Results**

We can first look at the Amplitude graph of Handel piece in Figure 1 and the frequency information in Figure 2. 

**Handel Example**

**spectrogram**

We then apply Gabor filter using Algorithm 1 to obtain the spectrogram graph as in Figure~\ref{fig:handelspec}. In Figure~\ref{fig:handelspec}, one Gabor filter is illustrated as the red curve in subplot 1, and the filtered result is in subplot 2.

![figure 1](https://github.com/EchoRLiu/Gabor-Transform/blob/master/handelspec.jpg)
![figure 2](https://github.com/EchoRLiu/Gabor-Transform/blob/master/amplitude1.jpg)
![figure 3](https://github.com/EchoRLiu/Gabor-Transform/blob/master/handelspec2.jpg)

**Window Width**

We then start to explore how window width influences our result, or how the trade off works. In Figure~\ref{fig:handelspec2}, compared to Figure~\ref{fig:handelspec}, we can clearly see there exist less low frequency and the information is more clear about which frequency happens during which time period. While in Figure~\ref{fig:handelspec3}, we obtain much more low frequency information, but our time information is much reduced. This is a very good example of the trade off between time and frequency information.

![figure 3](https://github.com/EchoRLiu/Gabor-Transform/blob/master/handelspec3.jpg)
![figure 4](https://github.com/EchoRLiu/Gabor-Transform/blob/master/handelspec4.jpg)

**Delta t**

We then start to explore how $\Delta t$ influences our result. In Figure~\ref{fig:handelspec4}, compared to Figure~\ref{fig:handelspec}, we can see there is a very sharp lines between each time window, meaning there is large discontinuity between each window. This is the result of Under-sampling. In Figure~\ref{fig:handelspec5}, we obtain much more continuous information, but it's also worth pointing out the code will run much slower and the information might not add much value.

![figure 5](https://github.com/EchoRLiu/Gabor-Transform/blob/master/handelspec5.jpg)
![figure 6]()

**Different Gabor windows**

We first use Mexcian hat window with result in  Figure~\ref{fig:mexicanhat}, then change to Shannon window as in Figure~\ref{fig:shannon}, and finally we use Step function window as in Figure~\ref{fig:step}. Comparing Figure~\ref{fig:handelspec} with the three different results, we can see how the window function's shape changes the spectrograms. From step function to Gaussian, to Mexican hat, and to Shannon (pay attention to the different red curves and different filtered results in the subplots 2), more information is being consideration for each time step, and thus similar to the result of a larger window, more low frequency information is included but the time information is reduced.

![figure 7](https://github.com/EchoRLiu/Gabor-Transform/blob/master/mexicanhat.jpg)
![figure 8](https://github.com/EchoRLiu/Gabor-Transform/blob/master/shannon.jpg)
![figure 9](https://github.com/EchoRLiu/Gabor-Transform/blob/master/step.jpg)

**Piano and Recorder**

In this example, we can see the difference between Piano and Recorder, which is that recorder does not produce overtone, while piano does. We also see how Algorithm~\ref{alg:part2_1}
hopes clear the overtone in the piano version of the music piece and reproduce music notes.

The amplitude and frequency information of the two music pieces are in Figure~\ref{fig:piano1} and Figure~\ref{fig:recorder1}. The two spectrograms are presented in Figure~\ref{fig:pianospec} and Figure~\ref{fig:recorderspec}. We can see the frequency graph or the spectrogram of the recorder piece is cleaner and does not have obvious layers of frequencies appearing in spectrogram.

![figure 10](https://github.com/EchoRLiu/Gabor-Transform/blob/master/piano1.jpg)
![figure 11](https://github.com/EchoRLiu/Gabor-Transform/blob/master/recorder1.jpg)

![figure 12](https://github.com/EchoRLiu/Gabor-Transform/blob/master/pianospec.png)
![figure 13](https://github.com/EchoRLiu/Gabor-Transform/blob/master/recorderspec.png)


Figure~\ref{fig:Note_eg} gives an example of the cleaned frequency information at time $t = .5$. The frequency here is the center frequency at that time. According to the cleaned up frequency information, we can reconstruct the music note here as: B A G A BBB AAA BBB BAGA BBBB AA BAG and the $\Delta t$ is around .6.

![figure 14](https://github.com/EchoRLiu/Gabor-Transform/blob/master/Note_eg.jpg)

**Summary and Conclusions**

As we can see in this example, Gabor transform is very powerful in giving us both time and frequency information on non-stationary data. We can also notice the trade off, or rather, the Heisenberg Uncertainty here: the more certain we are about the time, the less frequency information we have; the more frequency information we have, the less certain we are about what time the frequency happens. A solution for this trade off is wavelet analysis, which gives us more information on both time and frequency.
