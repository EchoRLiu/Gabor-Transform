close all; clear all; clc

load handel
v = y'/2; % Transpose is necessary for applying to filter later.
v=v(1:length(v)-1); % We delete the end data to make sure the length is even.
figure(1)
plot((1:length(v))/Fs,v); % Plot the Amplitude.
xlabel('Time [sec]');
ylabel('Amplitude');
title('Signal of Interest, v(n)');

% p8 = audioplayer(v,Fs);
% playblocking(p8);

%%

L=length(v)/Fs; % We set L to be the seconds of how long the sound lasts.
n=length(v);
t2=linspace(0,L,n+1); 
t=t2(1:n); % We assume period here and only need to use the first n.
k=(2*pi/L)*[0:n/2-1 -n/2:-1]; % Construted to match with fft results.
% 2pi rescale is necessary here since we assume the period to be 2pi.
ks=fftshift(k); % Shift it back to the right mathematical result.

vt=fft(v); % vt should be shifted back if needed for plotting.

figure(2); 
subplot(2,1,1)
plot((1:n)/Fs, v, 'k'); 
% (1:n)/Fs would be the time, since we scaled n with Fs before.
set(gca,'Fontsize',[14]); xlabel('Time (t)'); ylabel('V(t)');

subplot(2,1,2)
plot(ks, abs(fftshift(vt))/max(abs(vt)),'k');
% abs is applied to avoid imaginary part, max is being divivded to rescale.
% fftshift is applied to make sure it's in the correct mathematical version.
set(gca,'Fontsize',[14]); xlabel('Frequency'); ylabel('FFT(V)');

%%

% Part 1.

% #1.
figure(3);
vgt_spec=[];
tslide=0:.1:8.9; % 8.9 is L, the total time of the piece.
for j=1:length(tslide)
    g=exp(-2*(t-tslide(j)).^2);
    vg=g.*v; % Piecewise multiplication.
    vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))]; % Add new row to vgt_spec.
    % fftshit is applied to make sure it's the right mathematical form.
    % abs is applied to avoid imaginary parts.
    subplot(3,1,1),plot(t,v,'k',t,g,'r');
    subplot(3,1,2),plot(t,vg,'k');
    % Look at the current frequencies during certain delta time.
    subplot(3,1,3),plot(ks,abs(fftshift(vgt))/max(abs(vgt)));
    % Use ks instead of k, 
    drawnow;
    pause(.1);
end

pcolor(tslide, ks, vgt_spec.'); 
% vgt_spec is 90by73112, we need to transpose it to correctly plot it.
% . is to make sure the starting time is correct.
shading interp;
set(gca,'Ylim', [-12000 12000],'Fontsize',[14]);
colormap(hot);

%%

% # 2.
figure(4);
vgt_spec=[];
tslide=0:.1:8.9;
for j=1:length(tslide)
    % We change the a to 50, which would be narrower window.
    g=exp(-50*(t-tslide(j)).^2);
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
    subplot(3,1,1),plot(t,v,'k',t,g,'r');
    subplot(3,1,2),plot(t,vg,'k');
    subplot(3,1,3),plot(ks,abs(fftshift(vgt))/max(abs(vgt)));
    drawnow;
    pause(.1);
end

pcolor(tslide, ks, vgt_spec.'); shading interp;
set(gca,'Ylim', [-12000 12000], 'Fontsize',[14]);
xlabel('Time'), ylabel('Frequency'), 
title('The spectrogram of Handel with b = 50');
colormap(hot);

% Less low frequency; more accurate time information.

%%

figure(5);
vgt_spec=[];
tslide=0:.1:8.9;
for j=1:length(tslide)
    % We change the a to .5.
    g=exp(-0.5*(t-tslide(j)).^2);
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
    subplot(3,1,1),plot(t,v,'k',t,g,'r');
    subplot(3,1,2),plot(t,vg,'k');
    subplot(3,1,3),plot(ks,abs(fftshift(vgt))/max(abs(vgt)));
    drawnow;
    pause(.1);
end

pcolor(tslide, ks, vgt_spec.'); shading interp;
set(gca,'Ylim', [-12000 12000],'Fontsize',[14]);
xlabel('Time'), ylabel('Frequency'), 
title('The spectrogram of Handel with b = .5');
colormap(hot);

%%

% # 3.
figure(6);
vgt_spec=[];
tslide=0:1.0:8.9; % Under-sampling.
for j=1:length(tslide)
    g=exp(-2.*(t-tslide(j)).^2);
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
    subplot(3,1,1),plot(t,v,'k',t,g,'r');
    subplot(3,1,2),plot(t,vg,'k');
    subplot(3,1,3),plot(ks,abs(fftshift(vgt))/max(abs(vgt)));
    drawnow;
    pause(.1);
end

pcolor(tslide, ks, vgt_spec.'); shading interp;
set(gca,'Ylim', [-12000 12000],'Fontsize',[14]);
xlabel('Time'), ylabel('Frequency'), 
title('The spectrogram of Handel with delta t = 1.0');
colormap(hot);

%%

figure(7);
vgt_spec=[];
tslide=0:.05:8.9; % Over-sampling.
for j=1:length(tslide)
    g=exp(-2.*(t-tslide(j)).^2);
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
    subplot(3,1,1),plot(t,v,'k',t,g,'r');
    subplot(3,1,2),plot(t,vg,'k');
    subplot(3,1,3),plot(ks,abs(fftshift(vgt))/max(abs(vgt)));
    drawnow;
    pause(.1);
end

pcolor(tslide, ks, vgt_spec.'); shading interp;
set(gca,'Ylim', [-12000 12000],'Fontsize',[14]);
xlabel('Time'), ylabel('Frequency'), 
title('The spectrogram of Handel with delta t = .05');
colormap(hot);

%%

% # 3.
% Start with Gaussian filter.
figure(8);
vgaut_spec=[];
tslide=0:.1:8.9;
for j=1:length(tslide)
    gaussian=exp(-20.0*(t-tslide(j)).^2);
    vgau=gaussian.*v; vgaut=fft(vgau);
    vgaut_spec=[vgaut_spec; abs(fftshift(vgaut))];
    subplot(3,1,1),plot(t,v,'k',t,gaussian,'r');
    subplot(3,1,2),plot(t,vgau,'k');
    subplot(3,1,3),plot(ks,abs(fftshift(vgaut))/max(abs(vgaut)));
    drawnow;
    pause(.1);
end

pcolor(tslide, ks, vgaut_spec.'); shading interp;
set(gca,'Ylim', [-12000 12000],'Fontsize',[14]);
xlabel('Time'), ylabel('Frequency'), 
title('The spectrogram of Handel with Gaussian filter');
colormap(hot);

% Signal concentrated in time must have a spread out Fourier Transform,
% correlating with wide range of frequencies; 

%%

% Mexican Hat filter.
figure(9)
vgt_spec=[];
tslide=0:.1:8.9;
sig = .3;% sqrt(1/2); .1

for j=1:length(tslide)
    g=(1-((t-tslide(j))/sig).^2).*exp(-(t-tslide(j)).^2/(2*sig^2))*2/(sqrt(3*sig)*pi^(1/4));
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
    subplot(3,1,1),plot(t,v,'k',t,g,'r');
    subplot(3,1,2),plot(t,vg,'k');
    subplot(3,1,3),plot(ks,abs(fftshift(vgt))/max(abs(vgt)));
    drawnow;
    pause(.1);
end

pcolor(tslide, ks, vgt_spec.'); shading interp;
set(gca,'Ylim', [-12000 12000],'Fontsize',[14]);
xlabel('Time'), ylabel('Frequency'), 
title('The spectrogram of Handel with Mexican hat filter');
colormap(hot);

%%

% Shannon Wavelet filter.
figure(10)
vgt_spec=[];
tslide=0:.1:8.9;

for j=1:length(tslide)
    g=sinc((t-tslide(j))/2).*cos(3*pi*(t-tslide(j))/2);
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
    subplot(3,1,1),plot(t,v,'k',t,g,'r');
    subplot(3,1,2),plot(t,vg,'k');
    subplot(3,1,3),plot(ks,abs(fftshift(vgt))/max(abs(vgt)));
    drawnow;
    pause(.1);
end

pcolor(tslide, ks, vgt_spec.'); shading interp;
set(gca,'Ylim', [-12000 12000],'Fontsize',[14]);
xlabel('Time'), ylabel('Frequency'), 
title('The spectrogram of Handel with Shannon filter');
colormap(hot);

%%

% Step function.
figure(11)
vgt_spec=[];
slide=.1; % In case we want to change the slide.
tslide=0:slide:8.9;

for j=1:length(tslide)
    g=step_func(slide, tslide(j),t,n,L); % Please see the function definded at the end.
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
    subplot(3,1,1),plot(t,v,'k',t,g,'r');
    subplot(3,1,2),plot(t,vg,'k');
    subplot(3,1,3),plot(ks,abs(fftshift(vgt))/max(abs(vgt)));
    drawnow;
    pause(.1);
end

pcolor(tslide, ks, vgt_spec.'); shading interp;
set(gca,'Ylim', [-12000 12000],'Fontsize',[14]);
xlabel('Time'), ylabel('Frequency'), 
title('The spectrogram of Handel with step function filter');
colormap(hot);

%%

% Part 2.

close all; clear all; clc

tr_piano=16; % record time in seconds
y=audioread('/Users/yuhongliu/Downloads/music1.wav'); 
Fs=length(y)/tr_piano;

figure(1)
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)'); drawnow

%%
p8 = audioplayer(y,Fs); playblocking(p8);

%%

% Gabor filter for music1.wav

v = y'/2; % length of v is even.
L=length(v)/Fs; % We set L to be the seconds of how long the sound lasts.
n=length(v);
t2=linspace(0,L,n+1); 
t=t2(1:n);
k=(2*pi/L)*[0:n/2-1 -n/2:-1]; 
ks=fftshift(k);

vt=fft(v);

figure(2); subplot(2,1,1);
plot((1:n)/Fs, v, 'k');
set(gca,'Fontsize',[14]); xlabel('Time (t)'); ylabel('V(t)');

subplot(2,1,2);
plot(ks, abs(fftshift(vt))/max(abs(vt)),'k');
set(gca,'Fontsize',[14]); xlabel('Frequency'); ylabel('FFT(V)');

%%

vgt_spec=[];
tslide=0:.1:(length(v)/Fs);
omegas = [];
width = [.1 .5 5. 10. 50.];
for i = 1
    for j=1:length(tslide)
        g=exp(-width(i)*(t-tslide(j)).^2);
        vg=g.*v; vgt=fft(vg);
        abs_s_vgt = abs(fftshift(vgt));
        vgt_spec=[vgt_spec; abs_s_vgt];
        [pos_t, pos_ks] = ind2sub(size(abs_s_vgt), find(abs_s_vgt == max(abs_s_vgt)));
        omegas(j,1)=ks(1, pos_ks(2));
        % Since FFT returns both positive and negative frequencies, we only
        % need the positiv part here, hence the (2).
        filter_f = exp(-.2*(k-omegas(j,1)).^2);
        vgft = filter_f.*vgt;
        % vgf = ifft(vgft);
        figure(i+2)
        plot(ks, abs(fftshift(vgft))/max(abs(vgft)))
        set(gca, 'Xlim',[1400 2200]),
        xlabel('Frequency'),
        title(['The time is: ' num2str(j*.1)]),
        pause(.1);

        %subplot(3,1,1),plot(t,v,'k',t,g,'r');
        %subplot(3,1,2),plot(t,vg,'k');
        %subplot(3,1,3),plot(ks,abs(fftshift(vgt))/max(abs(vgt)));
        %drawnow;
    end
end

%%
pcolor(tslide, ks, vgt_spec.'); shading interp;
set(gca,'Ylim', [-60000 60000] , 'Fontsize',[14]);
colormap(hot);

%%

close all; clear all; clc

figure(1)
tr_rec=14; % record time in seconds
y=audioread('/Users/yuhongliu/Downloads/music2.wav'); 
Fs=length(y)/tr_rec;
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (recorder)');

% p8 = audioplayer(y,Fs); playblocking(p8);

%%

% Gabor filter for music2.wav

v = y'/2; % length of v is even.
L=tr_rec; % We set L to be the seconds of how long the sound lasts.
n=length(v);
t2=linspace(0,L,n+1); 
t=t2(1:n);
k=(2*pi/L)*[0:n/2-1 -n/2:-1]; 
ks=fftshift(k);

vt=fft(v);

figure(2); subplot(2,1,1);
plot((1:n)/Fs, v, 'k');
set(gca,'Fontsize',[14]); xlabel('Time (t)'); ylabel('V(t)');

subplot(2,1,2);
plot(ks, abs(fftshift(vt))/max(abs(vt)),'k');
set(gca,'Fontsize',[14]); xlabel('Frequency'); ylabel('FFT(V)');

%%

vgt_spec=[];
tslide=0:.1:(length(v)/Fs);
omegas = [];
width = [.1 .5 5.0 10. 50.];
for i = 1
    for j=1:length(tslide)
        g=exp(-width(i)*(t-tslide(j)).^2);
        vg=g.*v; vgt=fft(vg);
        abs_s_vgt = abs(fftshift(vgt));
        vgt_spec=[vgt_spec; abs_s_vgt];
        [pos_t, pos_ks] = ind2sub(size(abs_s_vgt), find(abs_s_vgt == max(abs_s_vgt)));
        omegas(j,1)=ks(1, pos_ks(2));
        % Since FFT returns both positive and negative frequencies, we only
        % need the positiv part here, hence the (2).
        filter_f = exp(-.2*(k-omegas(j,1)).^2);
        vgft = filter_f.*vgt;
        % vgf = ifft(vgft);
        figure(i+2)
        plot(ks, abs(fftshift(vgft))/max(abs(vgft)))
        set(gca, 'Xlim',[4000 8000]),
        xlabel('Frequency'),
        title(['The time is: ' num2str(j*.1)]),
        %pause(.1);

        %subplot(3,1,1),plot(t,v,'k',t,g,'r');
        %subplot(3,1,2),plot(t,vg,'k');
        %subplot(3,1,3),plot(ks,abs(fftshift(vgt))/max(abs(vgt)));
        %drawnow;
    end
end

%%

pcolor(tslide, ks, vgt_spec.'); shading interp;
set(gca, 'Fontsize',[14]);
colormap(hot);

%%

% Step function.
function step_g=step_func(slide, tslide, t, n, L)
[a,b] = size(t);
step_g = zeros(1,b);
if (t(end)-slide)>=tslide
    step_g((round(tslide*n/L)+1):round((tslide+slide)*n/L))=1.;
else
    step_g((round(tslide*n/L)+1):end)=1.;
end
end

