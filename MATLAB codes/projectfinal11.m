% MATLAB codes for comparing the performance of SINR vs SNR
clear all;

N = 20; %number of snapshots
M = 10; %number of sensors


tet_I= [30, 50]; % impinging angles of interfering sources

J  = length(tet_I); % Number of interfering sources

a_I = exp(1i*2.0*[0:M-1]'*pi*0.5*sin(tet_I*pi/180.0)); % interfering steering matrix
INR  = 1000; % interference noise ratio
SNR_dB = [-20 -19 -18 -17 -16 -15 -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5];

for i = 1:length(SNR_dB)
    SNR(i) = 10 ^ (SNR_dB(i)/10);
end

tet_no = 3; % nominal impinging angle
tet_ac = 5; % actual impinging angle(In Example 2)
R_I_n=INR*a_I*a_I'+eye(M,M); % ideal interference-plus-noise covariance matrix for optimal SINR

SINR_SMI=zeros(1,length(SNR));
SINR_LSMI=zeros(1,length(SNR));


for loop = 1:10
     a = exp(1i*2.0*[0:M-1]'*pi*0.5*sin(tet_no*pi/180.0)); %generation of actual signal steering vector a_tilde (Example 1)
     a_t = a; %Example 1
     %a_t = exp(1i*2.0*[0:M-1]'*pi*0.5*sin(tet_ac*pi/180.0)); %Example 2
     a_t = sqrt(10)*a_t/norm(a_t)';
     Ncount=0;
     for ps=SNR
         Ncount = Ncount + 1;
         R_h = zeros(M);
         for i = 1:N
          s = (randn(1,1) + 1i*randn(1,1))*sqrt(ps)*sqrt(0.5);  % signal waveform
          Is = (randn(J,1) + 1i*randn(J,1)).*sqrt(INR)'*sqrt(0.5); % interference waveform
          n = (randn(M,1) + 1i*randn(M,1))*sqrt(0.5); % noise
          x = a_I*Is + s*a_t + n; % snapshot with signal
          R_h = R_h + x*x';
         end
         R_h = R_h/N; %norm of covariance matrix
         
         R_dl_h = R_h + eye(M,M)*10; %loaded sample matrix inversion (LSMI)
         w_LSMI = inv(R_dl_h)*a;
         
         w_SMI = inv(R_h)*a; %sample matrix inversion (SMI)
         
         R_h = R_h + eye(M,M)*0.01; % diagonal loading
         U=chol(R_h);  % Cholesky factorization
         U_h=[real(U),-imag(U);imag(U),real(U)];  % 'U_bowl' defined based on 'U'
         M2 = 2*M;
         a_h=[real(a);imag(a)];  % 'a_bowl' difined based on 'a'
         a_b = [imag(a);-real(a)];
         FT = [1,zeros(1,M2);zeros(M2,1),U_h;0,a_h';zeros(M2,1),3*eye(M2,M2);0,a_b'];
         cvx_begin
            variable y(M2+1,1)
            minimize (y(1))
            c = FT*y;
            subject to
              c(1) >= norm(c(2:M2+1));
              c(M2+2)-1 >= norm(c(M2+3:2*M2+2));
              c(2*M2+3) == 0;
         cvx_end
         SINR_SMI_v(Ncount) = ps * (abs(w_SMI'*a_t) * abs(w_SMI'*a_t)) / abs(w_SMI'*R_I_n*w_SMI);  % SINR
         SINR_LSMI_v(Ncount) = ps * (abs(w_LSMI'*a_t) * abs(w_LSMI'*a_t)) / abs(w_LSMI'*R_I_n*w_LSMI);
         SINR_opt(Ncount) = ps * a_t' * inv(R_I_n) * a_t;
     end
     SINR_SMI     = ((loop-1)/loop)*SINR_SMI + (1/loop)*SINR_SMI_v; % just average
     SINR_LSMI    = ((loop-1)/loop)*SINR_LSMI + (1/loop)*SINR_LSMI_v;
end
save('Fig1_for_EX1_final','SINR_opt','SINR_SMI', 'SINR_LSMI') % Ex1
%save('Fig1_for_EX2_final','SINR_opt','SINR_SMI', 'SINR_LSMI') % Ex2

