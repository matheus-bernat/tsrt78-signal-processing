function [Phi,f]=sig2welch(x,varargin)
%Spectral analysis using Welch
[X,f]=batchdft(x,varargin{:});
Phi=mean(X,2);
