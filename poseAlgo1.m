close all; clear all;
%% Get files and do input checks
folderAll  = {'.\Artifact1'};
sensors = {'Sensor1'}
for ii = 1:length(folderAll)
    folder1 = folderAll{ii};

    files1 = dir(fullfile(folder1,'*.txt'));
    N = length(files1);
    if N == 0
        error('No files found in the folder');
    end

    %% Perform SVD
    for jj = 1:N
        fname1 = fullfile(folder1,files1(jj).name);
        dataT = dlmread(fname1);
        data1 = dataT(:,1:3); clear dataT; %Ignoring intensity/color data now.
        len1(jj) = length(data1);

        [r1,s1,v1] = svd(cov(data1));
        rAll(jj,:,:)=r1; %Store all the rotation matrices
        v2(jj,:) = v1(:,3); %Store all the plane normal vectors (dominant vector)
    end
    Rsum = squeeze(mean(rAll)); %Get the mean of the rotation matrix
    Rbar = getMeanRotFromSum(Rsum); % from Marek: this is correct mean rotation

    %% Perform checks on
    disp('Checking if rotation matrix is orthogonal')
    checkOk = Rbar*Rbar' - eye(3)
    vbar = mean(v2); %Get the mean of the plane normal vectors
    err1 = median(abs(checkOk(:)));
    if err1 < 1E-6
        disp('OK - Rotation matrix is orthogonal'),
    else
        error('Rotation matrix is NOT orthogonal'),
    end

    % Calculate the angle/alpha
    for jj = 1:length(files1)
        rm = squeeze(rAll(jj,:,:));
        deltaRM = Rbar*rm';
        alpha1(jj) = acos(0.5*trace(deltaRM)-0.5);
    end

    U = 95;
    P1 = prctile(alpha1,U);
    MULT = 1000;
    fprintf('Range of angles = %2.3f - %2.3f milliradian, %2.1f percentile = %2.3f milliradian\n',min(alpha1)*MULT,max(alpha1)*MULT, U, P1*MULT)

    %Plotting histogram
    figure(ii);clf
    hh = histogram(alpha1,20); hold on;
    y1 = 0:max(hh.Values);
    x1 = y1*0 + P1;

    plot(x1,y1,'r')
    legend('',sprintf('%dth percentile',U))
    title(sprintf('%s: Histogram and %2.1f percentile', sensors{ii},U))
    xlabel('Angle in radians'); ylabel('Frequency')
    saveas(gcf,strcat('MAT_',sensors{ii},'.png'))
    %break;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Rbar ] = getMeanRotFromSum( rSum )
%getMeanRotFromSum Ravg = getMeanRotFromSum( rSum )
%   sumR is 3x3 mat = simple algebraic sum / N
% Ravg is rotation matrix
%   the averaged calculated as in Eq.(3.7) Moakher 2002 "Means & Averaging
%   in the Group of Rotations', SIAM J.Matrix Analysis and Applications,
%   vol.24, pp.1-16

%Ravg = zeros(3,3);
RC = rSum'*rSum;
[U,D,V] = svd(RC);
dd = sqrt(ones(3,1) ./ diag(D));
Lambda = diag(dd);
Rbar = rSum*U*Lambda*U';
end



