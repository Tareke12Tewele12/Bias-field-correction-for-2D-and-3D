% This Matlab file demomstrates the method for simultaneous segmentation and bias field correction
% in Chunming Li et al's paper:
%    "Multiplicative intrinsic component optimization (MICO) for MRI bias field estimation and tissue segmentation",
%     Magnetic Resonance Imaging, vol. 32 (7), pp. 913-923, 2014
% Author: Chunming Li, all rights reserved
% E-mail: li_chunming@hotmail.com
% URL:  http://imagecomputing.org/~cmli/

clc;close all;clear all;
iterNum = 6;
N_region=2;  q=1.5;
%Img=imread('brainweb64.tif');
Nif =niftiread('D:\Spain_Udg\MiSa\Lab_\l1_preprocessing\braindata\t1_icbm_normal_1mm_pn0_rf20.nii');
yes =Nif(:,:,29)
save('D:\Spain_Udg\MiSa\Lab_\l1_preprocessing\source code\MICO_v0\MICO_v0\MICO_2D')
imshow(yes,[])
[rows, columns, numberOfColorChannels] = size(yes)

%yes=imread('D:\Spain_Udg\MiSa\Lab_\l1_preprocessing\source code\MICO_v0\MICO_v0\MICO_2D\brainweb59.tif');

%yes=imread('mprage171.tif');
yes = double(yes(:,:,1));
%load ROI
A=255;
yes_original = yes;
[nrow,ncol] = size(yes);n = nrow*ncol;

ROI = (yes>20); ROI = double(ROI);

tic

Bas=getBasisOrder3(nrow,ncol);
N_bas=size(Bas,3);
for ii=1:N_bas
    yesG{ii} = yes.*Bas(:,:,ii).*ROI;
    for jj=ii:N_bas
        GGT{ii,jj} = Bas(:,:,ii).*Bas(:,:,jj).*ROI;
        GGT{jj,ii} = GGT{ii,jj} ;
    end
end


energy_MICO = zeros(3,iterNum);

b=ones(size(yes));
for ini_num = 1:1
    C=rand(3,1);
    C=C*A;
    M=rand(nrow,ncol,3);
    a=sum(M,3);
    for k = 1 : N_region
        M(:,:,k)=M(:,:,k)./a;
    end
    
    [e_max,N_max] = max(M,[], 3);
    for kk=1:size(M,3)
        M(:,:,kk) = (N_max == kk);
    end
    
    M_old = M; chg=10000;
    energy_MICO(ini_num,1) = get_energy(yes,b,C,M,ROI,q);
    
    
    for n = 2:iterNum
        pause(0.1)
        
        [M, b, C]=  MICO(yes,q,ROI,M,C,b,Bas,GGT,yesG,1, 1);
        energy_MICO(ini_num,n) = get_energy(yes,b,C,M,ROI,q);
        
        figure(2),
        if(mod(n,1) == 0)
            PC=zeros(size(yes));
            for k = 1 : N_region
                PC=PC+C(k)*M(:,:,k);
            end
            subplot(241),imshow(uint8(yes)),title('original')
            subplot(242),imshow(PC.*ROI,[]); colormap(gray);
            iterNums=['segmentation: ',num2str(n), ' iterations'];
            title(iterNums);
            subplot(243),imshow(b.*ROI,[]),title('bias field')
            yes_bc = yes./b;  % bias field corrected image
            subplot(244),imshow(uint8(yes_bc.*ROI),[]),title('bias corrected')
            subplot(2,4,[5 6 7 8]),plot(energy_MICO(ini_num,:))
            xlabel('iteration number');
            ylabel('energy');
            pause(0.1)
        end
    end
end

[M,C]=sortMemC(M,C);
seg=zeros(size(yes));
for k = 1 : N_region
    seg=seg+k*M(:,:,k);   % label the k-th region 
end
figure;
subplot(141),imshow(yes,[]),title('Original image');
subplot(142),imshow(seg.*ROI,[]),title('Segmentation result');
subplot(143),imshow(b.*ROI,[]),title('bias field')
subplot(144),imshow(uint8(yes_bc.*ROI),[]),title('bias corrected')



