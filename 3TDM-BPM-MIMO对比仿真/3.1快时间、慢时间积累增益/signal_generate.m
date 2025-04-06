function [echo_baseband] = signal_generate(amp, range_round_trip, radar_parameter)
k = radar_parameter.K;
sample_num = radar_parameter.sample_num;
fc = radar_parameter.fc*ones(sample_num,1);
fast_time = ((0:sample_num-1)/radar_parameter.fs)';
delay = range_round_trip/radar_parameter.c;
echo_baseband = (amp * exp(1i*2*pi*(k*fast_time*delay + fc*delay  - 1/2*k*(delay).^2)));
% y = exp(1i*2*pi*fc*(tnrn+n*Ta-tobject)).*exp(1i*pi*k*(tnrn-tobject).^2) ;
% echo_baseband(:,i) = amp * exp(1i*2*pi*(k*fast_time*delay + delay * fc - 1/2*k*delay.^2));
end