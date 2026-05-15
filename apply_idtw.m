function [new_data] = apply_idtw(pre_data, fs)
% Input
%   - pre_data: three-dimensional data [number of channels, number of samples, number of trials]
%   - fs: sampling rate (Hz)
% Output
%   - new_data: data processed via the proposed IDTW [number of channels, number of samples, number of trials]
    new_data = zeros(size(pre_data));

    for i_ch=1:2
        flag_end = 0;
        count = 1;
        temp_sig = permute(squeeze(pre_data(i_ch,:,:)),[2 1]);
        pre_mean_sig = mean(temp_sig,1);

        while(flag_end == 0)
            [temp_sig2,~]= make_dtw2(temp_sig,fs,pre_mean_sig);
            current_mean_sig = mean(temp_sig2,1);

            if (sqrt(sum(current_mean_sig.^2)) * 0.01 > sqrt(sum((current_mean_sig - pre_mean_sig).^2))) || (count>=50)
                flag_end = 1;
                break;
            end

            pre_mean_sig = current_mean_sig;
            count = count + 1;
        end

        new_data(i_ch,:,:) = permute(temp_sig2,[2 1]);
    end
end

function [sig_r,indx] = make_dtw2(sig,fs,sig_m)
    sig_r = zeros(size(sig));
    n = length(sig_m);
    indx = zeros(size(sig));

    for i = 1:size(sig,1)
        c = zeros(n,n);
        sig_t = sig(i,:);

        for i1=1:n
            for i2=1:n
                c(i1,i2) = abs(sig_m(i1) - sig_t(i2));
            end
        end

        [dist,ix,iy] = dtw(sig_m,sig_t,'absolute');
        
        sig_r(i,1) = sig_t(1);
        indx(i,1) = 1;
        count = 1;
        for ii = 2:length(ix)
            d = ix(ii) - ix(ii-1);
            if d ~= 0
                count = count + 1;
                sig_r(i,count) = sig_t(iy(ii));
                indx(i,count) = iy(ii);
            end
        end

         n=3;
         wn=30;
         fn=fs/2;
         ftype='low';
         [b, a]=butter(n,wn/fn,ftype);
         sig_r(i,:)=filtfilt(b,a,sig_r(i,:));
    end
end