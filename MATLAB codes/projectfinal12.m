clear all;

M=10;  % No. of sensors
SNR_dB = -10;
ps = 10^(SNR_dB/10);

t_I= [30, 50]; % impinging angles of interfering sources
J  = length(t_I); % No. of interfering sources

INR  = 1000; % interference noise ratio


a_I = exp(1i*2.0*[0:M-1]'*pi*0.5*sin(t_I*pi/180.0)); % interfering steering matrix
R_I_n=INR*a_I*a_I'+eye(M,M); % ideal interference-plus-noise covariance matrix for computing optimal SINR

N=[1 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100];
SMI=zeros(1,length(N));
LSMI=zeros(1,length(N));

tet_no=3; % nominal impinging angle
tet_ac=5; % actual impinging angle(In Example 2)

for loop = 1:100;
    a = exp(1i*2.0*[0:M-1]'*pi*0.5*sin(tet_no*pi/180.0)); 
    %a_t = a;
    a_t = exp(1i*2.0*[0:M-1]'*pi*0.5*sin(tet_ac*pi/180.0)); %Example 2
    a_t = sqrt(10)*a_t/norm(a_t)';
    Ncount=0;
  for N=[1 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100];
      Ncount=Ncount+1;
      R_h = zeros(M);
      for l=1:N
          sw=(randn(1,1) + 1i*randn(1,1))*sqrt(ps)*sqrt(0.5);  % signal waveform
          Is = (randn(J,1) + 1i*randn(J,1)).*sqrt(INR)'*sqrt(0.5); % interference waveform
          n = (randn(M,1) + 1i*randn(M,1))*sqrt(0.5); % noise
          x = a_I*Is + sw*a_t + n; % snapshot with signal
          R_h = R_h + x*x';
      end
      R_h=R_h/length(N);
      
      w_SMI=inv(R_h)*a;
      
      R_dl_h = R_h + eye(M,M)*10;
      w_LSMI=inv(R_dl_h)*a;
      
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
      SMI_v(Ncount) = ps * (abs(w_SMI'*a_t) * abs(w_SMI'*a_t)) / abs(w_SMI'*R_I_n*w_SMI);  % SINR
      LSMI_v(Ncount) = ps * (abs(w_LSMI'*a_t) * abs(w_LSMI'*a_t)) / abs(w_LSMI'*R_I_n*w_LSMI);
      SINR_opt1(Ncount) = ps * a_t' * inv(R_I_n) * a_t;
  end
  SMI     = ((loop-1)/loop)*SMI + (1/loop)*SMI_v; % just average
  LSMI    = ((loop-1)/loop)*LSMI + (1/loop)*LSMI_v;
end
save('Fig2_for_EX2_final1','SINR_opt1','SMI', 'LSMI')
