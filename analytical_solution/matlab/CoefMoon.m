function [J,C,S,C_mtx,S_mtx] = CoefMoon
%COEFMOON Lunar gravity coefficients used by the analytical solution.
%
%   [J,C,S,C_mtx,S_mtx] = CoefMoon
%
%   Returns the zonal harmonics Jn (n=2..23) and, from data/Coef_moon.mat,
%   the fully denormalized tesseral/sectorial coefficients C_mtx, S_mtx
%   (degree 2..50, order 0..degree) used by RACook_const.m and
%   coeff_odd.m. C and S are a legacy subset (n=m only, degree 2..23) kept
%   for backward compatibility; C_mtx/S_mtx are what the solution actually
%   uses.
%
% normalised zonal coeff (from GRIM)
C20 = -0.908833998020E-04;
C22 = 0.346733143947E-04;
S22 = -0.392057200495E-10;
C30 = -0.319752790393E-05;
C33 = 0.122751621027E-04;
S33 = -0.177442570292E-05;
C40 = 0.323481934886E-05;
C44 = -0.600694789318E-05;
S44 = 0.392673085841E-05;
C50 = -0.223787885166E-06;
C55 = 0.311087474731E-05;
S55 = -0.275457108429E-05;
C60 = 0.381843306294E-05;
C66 = -0.468411670641E-05;
S66 = 0.722996725837E-05;
C70 = 0.559339280607E-05;
C77 = -0.182030463112E-05;
S77 = 0-.160008038436E-05;
C80 = 0.234683041324E-05;
C88 = -0.248561579818E-05;
S88 = 0.211650771319E-05;
C90 = -0.353091051156E-05;
C99 = -0.938457887352E-06;
S99 = 0.248814382461E-05;
C100 = -0.106930657929E-05;
C1010 = 0.947913756735E-06;
S1010 = -0.171933176281E-05;
C110 = -0.885206325693E-06;
C1111 = -0.288888330260E-05;
S1111 = -0.174634261598E-05;
C120 = -0.201225870859E-05; 
C1212 = 0.302810668694E-06;
S1212 = 0.124695798353E-05;
C130 =  0.256479133203E-06;
C1313 = 0.249932734075E-05;
S1313 = -0.266878879331E-05;
C140 = 0.492249665248E-06;
C1414 = -0.613492433245E-06;
S1414 = 0.107106362472E-05;
C150 = -0.228086318837E-06;
C1515 = 0.548312371917E-06;
S1515 = 0.777863375286E-06;
C160 = 0.556851031937E-06;
C1616 = -0.905389352854E-06;
S1616 = -0.918816470920E-06;
C170 = -0.103303359350E-05;
C1717 = -0.380714295482E-06;
S1717 = 0.150465577167E-05;
C180 = -.520320699644E-06;
C1818 = -.209680416757E-06;
S1818 = -0.285137994367E-06;
C190 = 0.215655284627E-06;
C1919 = 0.602056696503E-06;
S1919 = 0.365326151124E-06;
C200 = 0.335562704824E-06;
C2020 = -0.135259305585E-06;
S2020 = -0.358791137249E-06;
C210 = -0.560127670919E-07;
C2121 = -0.334111629823E-06;
S2121 = 0.324093861714E-06;
C220 = -0.841390472495E-07;
C2222 = 0.356170428935E-06;
S2222 = 0.656489854634E-06;
C230 = -0.385438615758E-06;
C2323 = -0.106136625412E-05;
S2323 = 0.249418624727E-06;
J2=-C20.*sqrt(1+2.*2);
J3=-C30.*sqrt(1+3.*2);
J4=-C40.*sqrt(1+4.*2);
J5=-C50.*sqrt(1+5.*2);
J6=-C60.*sqrt(1+6.*2);
J7=-C70.*sqrt(1+7.*2);
J8=-C80.*sqrt(1+8.*2);
J9=-C90.*sqrt(1+9.*2);
J10=-C100.*sqrt(1+10.*2);
J11=-C110.*sqrt(1+11.*2);
J12=-C120.*sqrt(1+12.*2);
J13=-C130.*sqrt(1+13.*2);
J14=-C140.*sqrt(1+14.*2);
J15=-C150.*sqrt(1+15.*2);
J16=-C160.*sqrt(1+16.*2);
J17=-C170.*sqrt(1+17.*2);
J18=-C180.*sqrt(1+18.*2);
J19=-C190.*sqrt(1+19.*2);
J20=-C200.*sqrt(1+20.*2);
J21=-C210.*sqrt(1+21.*2);
J22=-C220.*sqrt(1+22.*2);
J23=-C230.*sqrt(1+23.*2);

%factorial(n+m)/(factorial(n-m)*(2*n+1)*2)
C22 =  C22/sqrt(factorial(4)/2/(2*2+1));
S22 = S22/sqrt(factorial(4)/2/(2*2+1));
C33 =  C33/sqrt(factorial(6)/2/(2*3+1));
S33 = S33/sqrt(factorial(6)/2/(2*3+1));
C44 = C44/sqrt(factorial(8)/2/(2*4+1));
S44 = S44/sqrt(factorial(8)/2/(2*4+1));
C55 = C55/sqrt(factorial(10)/2/(2*5+1));
S55 = S55/sqrt(factorial(10)/2/(2*5+1));
C66 =  C66/sqrt(factorial(12)/2/(2*6+1));
S66 = S66/sqrt(factorial(12)/2/(2*6+1));
C77 =  C77/sqrt(factorial(14)/2/(2*7+1));
S77 = S77/sqrt(factorial(14)/2/(2*7+1));
C88 = C88/sqrt(factorial(16)/2/(2*8+1));
S88 = S88/sqrt(factorial(16)/2/(2*8+1));
C99 = C99/sqrt(factorial(18)/2/(2*9+1));
S99 = S99/sqrt(factorial(18)/2/(2*9+1));
C1010 = C1010/sqrt(factorial(20)/2/(2*10+1));
S1010 = S1010/sqrt(factorial(20)/2/(2*10+1));
C1111 = C1111/sqrt(factorial(22)/2/(2*11+1));
S1111 = S1111/sqrt(factorial(22)/2/(2*11+1));
C1212 = C1212/sqrt(factorial(24)/2/(2*12+1));
S1212 = S1212/sqrt(factorial(24)/2/(2*12+1));
C1313 = C1313/sqrt(factorial(26)/2/(2*13+1));
S1313 = S1313/sqrt(factorial(26)/2/(2*13+1));
C1414 = C1414/sqrt(factorial(28)/2/(2*14+1));
S1414 = S1414/sqrt(factorial(28)/2/(2*14+1));
C1515 = C1515/sqrt(factorial(30)/2/(2*15+1));
S1515 = S1515/sqrt(factorial(30)/2/(2*15+1));
C1616 = C1616/sqrt(factorial(32)/2/(2*16+1));
S1616 = S1616/sqrt(factorial(32)/2/(2*16+1));
C1717 = C1717/sqrt(factorial(34)/2/(2*17+1));
S1717 = S1717/sqrt(factorial(34)/2/(2*17+1));
C1818 = C1818/sqrt(factorial(36)/2/(2*18+1));
S1818 = S1818/sqrt(factorial(36)/2/(2*18+1));
C1919 = C1919/sqrt(factorial(38)/2/(2*19+1));
S1919 = S1919/sqrt(factorial(38)/2/(2*19+1));
C2020 = C2020/sqrt(factorial(40)/2/(2*20+1));
S2020 = S2020/sqrt(factorial(40)/2/(2*20+1));
C2121 = C2121/sqrt(factorial(42)/2/(2*21+1));
S2121 = S2121/sqrt(factorial(42)/2/(2*21+1));
C2222 = C2222/sqrt(factorial(44)/2/(2*22+1));
S2222 = S2222/sqrt(factorial(44)/2/(2*22+1));
C2323 = C2323/sqrt(factorial(46)/2/(2*23+1));
S2323 = S2323/sqrt(factorial(46)/2/(2*23+1));
C = [0,C22,C33,C44,C55,C66,C77,C88,C99,C1010,C1111,C1212,C1313,C1414,C1515,C1616,C1717,C1818,C1919,C2020,C2121,C2222,C2323];
S = [0,S22,S33,S44,S55,S66,S77,S88,S99,S1010,S1111,S1212,S1313,S1414,S1515,S1616,S1717,S1818,S1919,S2020,S2121,S2222,S2323];
J = [0,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,J20,J21,J22,J23];

this_dir = fileparts(mfilename('fullpath'));
coef_candidates = {
    fullfile(this_dir, 'Coef_moon.mat'), ...
    fullfile(this_dir, 'Coef_moon'), ...
    fullfile(this_dir, '..', 'data', 'Coef_moon.mat'), ...
    fullfile(this_dir, '..', 'data', 'Coef_moon') ...
    };

coef_path = '';
for idx = 1:numel(coef_candidates)
    if exist(coef_candidates{idx}, 'file')
        coef_path = coef_candidates{idx};
        break
    end
end

if isempty(coef_path)
    error('CoefMoon:MissingFile', 'Unable to find Coef_moon data file.');
end

coef = load(coef_path);
C_mtx = zeros(49,51);
S_mtx = zeros(49,51);
for i = 2:50
    for j = 0:i
        if j == 0
            k = 1;
        else
            k = 2;
        end
        C_mtx(i-1,j+1) = coef.Clm(i-1,j+1)./sqrt(factorial(i+j)/(factorial(i-j)*(2*i+1)*k));
        S_mtx(i-1,j+1) = coef.Slm(i-1,j+1)./sqrt(factorial(i+j)/(factorial(i-j)*(2*i+1)*k));
    end
end

J = [0; -C_mtx(:,1)];
end