%% run frequency check
%  [status, dat, fitme, esf, nbin, del2] = sfrmat3(io, del, weight, a, oecfname);
%       From a selected edge area of an image, the program computes
%       the ISO slanted edge SFR. Input file can be single or
%       three-record file. Many image formats are supported. The image
%       is displayed and a region of interest (ROI) can be chosen, or
%       the entire field will be selected by not moving the mouse
%       when defining an ROI (simple click). Either a vertical or horizontal
%       edge features can be analized.
%  Input arguments:
%      io  (optional)
%        0 = (default) R,G,B,Lum SFRs + edge location(s)
%          = 'sfrmat2'  R,G,B,Lum SFRs + edge location(s)but
%            with the same calculations as the previous version, sfrmat2
%        1 = Non GUI usage with supplied data array
%      del (optional) sampling interval in mm or pixels/inch
%          If dx < 1 it is assumed to be sampling pitch in mm
%          If io = 1 (see below, no GUI) and del is not specified,
%          it is set equal to 1, so frequency is given in cy/pixel.
%      weight (optiona) default 1 x 3 r,g,b weighs for luminance weighting
%      a   (required if io =1) an nxm or nxmx3 array of data
%      oename  optional name of oecf LUT file containing 3xn or 1xn array
%
% Returns: 
%       status = 0 if normal execution
%       dat = computed sfr data
%       fitme = coefficients for the linear equations for the fit to
%               edge locations for each color-record. For a 3-record
%               data file, fitme is a (4 x 3) array, with the last column
%               being the color misregistration value (with green as 
%               reference).
%       esf = supersampled edge-spread functin array
%       nbin = binning factor used
%       del2 = sampling interval for esf, from which the SFR spatial
%              frequency sampling is was computed. This will be 
%              approximately  4  times the original image sampling.

  
%% run
cd('D:\Thesis_ light field_D\my works\Light_F_framework\Toolbox\Benchmarking\sfrmat3_post')
a=sfrmat3

%% plot
plot(freq,interpolated_Barycentric,'k',...
    freq,interpolated_bicubic,'r',...
    freq,interpolated_Gaussian,'g',...
    freq,interpolated_Linear,'b',...
    freq,interpolated_NN,'c',...
    freq,interpolated_Sinc,'m')
legend('Barycentric','Bicubic','Gaussian','Linear','NN','Sinc'); 