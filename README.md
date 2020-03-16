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
\subsubsection{Piano and Recorder}
"Mary had a little lamb" would be played on both piano and recorder. To see the difference between piano and recorder, two spectrograms are produced. We can see that compared to piano, recorder does not have the overtone structure, meaning there is no layers of frequencies appearing in the spectrogram. To obtain a good music score for the piano piece, GT will first be used to obtain the spectrogram for us to inspect which frequencies are producing overtones, then Fourier transform is used to filter out frequencies overtones.

%  Theoretical Background
\section{Theoretical Background}
As we learned from our textbook \cite{kutz_2013}, to solve the loss of time information problem, Gabor introduced the way to slice up time, use Gabor filter around certain time window, and then apply Fourier transform. Gabor transform, or short-time Fourier Transform is defined as follows:
\begin{equation}
    G[f](t,\omega) = \Tilde{f}_g(t, \omega)=\int^{\infty}_{-\infty}f(\tau)g(\tau-t)e^{-i\omega\tau}d\tau=(f,\bar{g}_{t,\omega})
    \label{eqn:fourierseries}
\end{equation}
This would be our main tool for this whole project. A main Gabor filter we would be using here is Gaussian window:
\begin{equation}
    g(t)=e^{-a(t-b)^2}
    \label{eqn:gaussianfilter}
\end{equation}
where a is the window width and b is the time point we can change to slide through the whole time duration, which we can adjust to extact time and frequency information.

To solve the drawback of GT, wavelet is introduced. The basic idea of wavelet is that we don't have to decide which a, or window width, we want to use when filtering. When using wavelet, a certain range of the window width would be used to extract more accurate information, thus resolving the trade off issue we would have when using Garbor transform. As illustrated in our textbook \cite{kutz_2013}, one good example is Mexican Hat wavelet:
\begin{equation}
    \psi(t)=(1-t^2)e^{-t^2/2}=-d^2/dt^2(e^{-t^2/2})
    \label{eqn:mexicanhat1}
\end{equation}
\begin{equation}
    \hat{\psi}(\omega)=\hat{\psi}_{1,0}(\omega)=\sqrt{2\pi}\omega^2e^{-\omega^2/2}
    \label{eqn:mexicanhat2}
\end{equation}

% Algorithm Implementation and Development
\section{Algorithm Implementation and Development}
Here, the way the code is developed is illustrated.
\begin{enumerate}
    \item handel music is loaded from Matlab to y, we then transpose it as store in v, which is a necessary step for applying Gabor filter later. The final point was taken out since we want to make sure the length is even for faster fft.
    \item the music time length, time domain and frequency domain is created.
    \item The following illustrate how Gabor filter is applied to the music piece to obtain spectrogram.
    
    \begin{algorithm}
    \begin{algorithmic}
    \STATE{Define \texttt{tslide} as a time vector from 0 to the time length of the music piece with distance $\Delta t$}
    \STATE{Define $vgt spec$ to store the Fourier transformed and filtered result.}
    \FOR{$j = 1:length(tslide)$}
        \STATE{Define $g$ as the Gaussian window for filtering out information outside of certain window width}
        \STATE{Multiply $g$ with v and store the result as $vg$}
        \STATE{Fourier Transform $vg$, $vgt = fft(vg)$}
        \STATE{Add the result to $vgt spec$.}
    \ENDFOR
    \STATE{Visualize $vgt spec$, which would give us the spectrogram.}
    \end{algorithmic}
    \caption{Gabor Transform}
    \label{alg:part1_1}
    \end{algorithm}
    
    \item The window width is then changed to $50$ and $.5$ in Algorithm~\ref{alg:part1_1} to see the difference of the spectrograms and analyze the trade off.
    \item The distance $\Delta t$ is then changed to $1.0$ and $.05$ in Algorithm~\ref{alg:part1_1} to see the difference of under-sampling and over-sampling.
    \item The filter function is then changed to illustrate different Gabor windows in Algorithm~\ref{alg:part1_1}, such as Gaussian window as in Algorithm~\ref{alg:part1_1}, Mexican hat wavelet, Shannon wavelet, step function. Different filter function $g$ is provided here.
    \begin{equation}
    g_1 = e^{-b(t-tslide(j)^2)}
    \label{eqn:filter1}
    \end{equation}
    \begin{equation}
    g_2=(1-(\frac{t-tslide(j)}{b})^2)\frac{2}{\sqrt{3b}\pi^{\frac{1}{4}}}e^{-\frac{(t-tslide(j))^2}{2b^2}}
    \label{eqn:filter2}
    \end{equation}
    \begin{equation}
    g_3=sinc(\frac{t-tslide(j)}{2})cos(3\pi\frac{t-tslide(j)}{2})
    \label{eqn:filter3}
    \end{equation}
    \begin{equation}
    g_4= \begin{cases} 
      1 & tslide(j) - \Delta t \leq t\leq tslide(j) \\
      0 & otherwise \\
    \end{cases}
    \label{eqn:filter4}
    \end{equation}
    \item The following would be for analyzing the piano and recorder cases. We first use Algorithm~\ref{alg:part1_1} above to create spectrograms for both pieces and observe their differences.
    \item Then, Gaussian filter is applied after Fourier Transform to filter out overtones using Algorithm~\ref{alg:part2_1}. We can get a very clean result showing us the music notes played during certain time period.
    
     \begin{algorithm}
    \begin{algorithmic}
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
    \end{algorithmic}
    \caption{Overtone Filter}
    \label{alg:part2_1}
    \end{algorithm}
    
\end{enumerate}

% Computational Results
\section{Computational Results}
We can first look at the Amplitude graph of Handel piece in Figure~\ref{fig:amplitude1} and the frequency information in Figure~\ref{fig:simplefft1}. 

\subsection{Handel Example}
\subsubsection{spectrogram}
We then apply Gabor filter using Algorithm~\ref{fig:part1_1} to obtain the spectrogram graph as in Figure~\ref{fig:handelspec}. In Figure~\ref{fig:handelspec}, one Gabor filter is illustrated as the red curve in subplot 1, and the filtered result is in subplot 2.

\begin{figure}
    \centering
    \includegraphics[width=0.5\linewidth]{amplitude1.jpg}
    \caption{Amplitude of Handel piece}
    \label{fig:amplitude1}
\end{figure}
\begin{figure}
    \centering
    \includegraphics[width=0.5\linewidth]{simplefft1.jpg}
    \caption{Frequency information of Handel}
    \label{fig:simplefft1}
\end{figure}
\begin{figure}
    \centering
    \includegraphics[width=0.5\linewidth]{handelspec.jpg}
    \caption{Spectrogram of Handel piece}
    \label{fig:handelspec}
\end{figure}

\subsubsection{Window Width}
We then start to explore how window width influences our result, or how the trade off works. In Figure~\ref{fig:handelspec2}, compared to Figure~\ref{fig:handelspec}, we can clearly see there exist less low frequency and the information is more clear about which frequency happens during which time period. While in Figure~\ref{fig:handelspec3}, we obtain much more low frequency information, but our time information is much reduced. This is a very good example of the trade off between time and frequency information.

\begin{figure}
    \centering
    \includegraphics[width=0.5\linewidth]{handelspec2.jpg}
    \caption{Spectrogram of Handel piece with small window width}
    \label{fig:handelspec2}
\end{figure}
\begin{figure}
    \centering
    \includegraphics[width=0.5\linewidth]{handelspec3.jpg}
    \caption{Spectrogram of Handel piece with large window width}
    \label{fig:handelspec3}
\end{figure}

\subsubsection{Delta t}
We then start to explore how $\Delta t$ influences our result. In Figure~\ref{fig:handelspec4}, compared to Figure~\ref{fig:handelspec}, we can see there is a very sharp lines between each time window, meaning there is large discontinuity between each window. This is the result of Under-sampling. In Figure~\ref{fig:handelspec5}, we obtain much more continuous information, but it's also worth pointing out the code will run much slower and the information might not add much value.

\begin{figure}
    \centering
    \includegraphics[width=0.5\linewidth]{handelspec4.jpg}
    \caption{Spectrogram of Handel piece with large delta t - Under-sampling}
    \label{fig:handelspec4}
\end{figure}
\begin{figure}
    \centering
    \includegraphics[width=0.5\linewidth]{handelspec5.jpg}
    \caption{Spectrogram of Handel piece with small delta t - Over-sampling}
    \label{fig:handelspec5}
\end{figure}

\subsubsection{Different Gabor windows}
We first use Mexcian hat window with result in  Figure~\ref{fig:mexicanhat}, then change to Shannon window as in Figure~\ref{fig:shannon}, and finally we use Step function window as in Figure~\ref{fig:step}. Comparing Figure~\ref{fig:handelspec} with the three different results, we can see how the window function's shape changes the spectrograms. From step function to Gaussian, to Mexican hat, and to Shannon (pay attention to the different red curves and different filtered results in the subplots 2), more information is being consideration for each time step, and thus similar to the result of a larger window, more low frequency information is included but the time information is reduced.

\begin{figure}
    \centering
    \includegraphics[width=0.5\linewidth]{mexicanhat.jpg}
    \caption{Spectrogram of Handel piece with Mexican Hat window}
    \label{fig:mexicanhat}
\end{figure}
\begin{figure}
    \centering
    \includegraphics[width=0.5\linewidth]{shannon.jpg}
    \caption{Spectrogram of Handel piece with Shannon window}
    \label{fig:shannon}
\end{figure}
\begin{figure}
    \centering
    \includegraphics[width=0.5\linewidth]{step.jpg}
    \caption{Spectrogram of Handel piece with Step function window}
    \label{fig:step}
\end{figure}

\subsection{Piano and Recorder}
In this example, we can see the difference between Piano and Recorder, which is that recorder does not produce overtone, while piano does. We also see how Algorithm~\ref{alg:part2_1}
hopes clear the overtone in the piano version of the music piece and reproduce music notes.

The amplitude and frequency information of the two music pieces are in Figure~\ref{fig:piano1} and Figure~\ref{fig:recorder1}. The two spectrograms are presented in Figure~\ref{fig:pianospec} and Figure~\ref{fig:recorderspec}. We can see the frequency graph or the spectrogram of the recorder piece is cleaner and does not have obvious layers of frequencies appearing in spectrogram.

\begin{figure}
    \centering
    \includegraphics[width=0.5\linewidth]{piano1.jpg}
    \caption{Amplitude and Frequency information for Piano piece}
    \label{fig:piano1}
\end{figure}
\begin{figure}
    \centering
    \includegraphics[width=0.5\linewidth]{recorder1.jpg}
    \caption{Amplitude and Frequency information for Recorder piece}
    \label{fig:recorder1}
\end{figure}

\begin{figure}
    \centering
    \includegraphics[width=0.5\linewidth]{pianospec.png}
    \caption{Spectrogram for Piano piece}
    \label{fig:pianospec}
\end{figure}
\begin{figure}
    \centering
    \includegraphics[width=0.5\linewidth]{recorderspec.png}
    \caption{Spectrogram for Recorder piece}
    \label{fig:recorderspec}
\end{figure}


Figure~\ref{fig:Note_eg} gives an example of the cleaned frequency information at time $t = .5$. The frequency here is the center frequency at that time. According to the cleaned up frequency information, we can reconstruct the music note here as: B A G A BBB AAA BBB BAGA BBBB AA BAG and the $\Delta t$ is around .6.

\begin{figure}
    \centering
    \includegraphics[width=0.5\linewidth]{Note_eg.jpg}
    \caption{Cleaned up frequency for piano piece}
    \label{fig:Note_eg}
\end{figure}

**Summary and Conclusions**

As we can see in this example, Gabor transform is very powerful in giving us both time and frequency information on non-stationary data. We can also notice the trade off, or rather, the Heisenberg Uncertainty here: the more certain we are about the time, the less frequency information we have; the more frequency information we have, the less certain we are about what time the frequency happens. A solution for this trade off is wavelet analysis, which gives us more information on both time and frequency.
