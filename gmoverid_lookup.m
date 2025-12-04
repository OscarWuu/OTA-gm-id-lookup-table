clear all, close all, clc

function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_60n(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L60n.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L60n_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end

function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_80n(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L80n.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L80n_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end

function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_100n(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L100n.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L100n_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end

function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_150n(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L150n.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L150n_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end

function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_200n(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L200n.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L200n_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end
function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_250n(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L250n.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L250n_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end
function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_300n(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L300n.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L300n_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end
function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_350n(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L350n.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L350n_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end
function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_400n(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L400n.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L400n_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end
function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_500n(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L500n.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L500n_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end
function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_600n(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L600n.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L600n_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end

function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_700n(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L700n.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L700n_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end
function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_800n(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L800n.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L800n_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end
function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_900n(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L900n.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L900n_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end
function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_1u(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L1u.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L1u_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end
function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_1d1u(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L1d1u.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L1d1u_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end
function [gmoverid_n, idw_n, gmdoverid_p, idw_p] = lookup_idw_1d2u(target_x) % L=200n
    % --- Load CSV file ---
    data1 = readmatrix('idw_vs_gmid_L1d2u.csv');   % change file name if needed
    data2 = readmatrix('idw_vs_gmid_L1d2u_p.csv'); 

    xn = data1(:,1);   % first column
    yn = data1(:,2);   % second column
    xp = data2(:,1);
    yp = data2(:,2);

    % --- Find the index of x closest to target_x ---
    [~, idxn] = min(abs(xn - target_x));
    [~, idxp] = min(abs(xp - target_x));

    % --- Return values ---
    gmoverid_n = xn(idxn);
    idw_n = yn(idxn);
    gmdoverid_p = xp(idxp);
    idw_p = abs(yp(idxp));
end


[a(1), b(1), c(1), d(1)] = lookup_idw_200n(1);
[a(2), b(2), c(2), d(2)] = lookup_idw_200n(2);
[a(3), b(3), c(3), d(3)] = lookup_idw_200n(3);
[a(4), b(4), c(4), d(4)] = lookup_idw_200n(4);
[a(5), b(5), c(5), d(5)] = lookup_idw_200n(5);
[a(6), b(6), c(6), d(6)] = lookup_idw_200n(6);
[a(7), b(7), c(7), d(7)] = lookup_idw_200n(7);
[a(8), b(8), c(8), d(8)] = lookup_idw_200n(8);
[a(9), b(9), c(9), d(9)] = lookup_idw_200n(9);
[a(10), b(10), c(10), d(10)] = lookup_idw_200n(10);
[a(11), b(11), c(11), d(11)] = lookup_idw_200n(11);
[a(12), b(12), c(12), d(12)] = lookup_idw_200n(12);
[a(13), b(13), c(13), d(13)] = lookup_idw_200n(13);
[a(14), b(14), c(14), d(14)] = lookup_idw_200n(14);
[a(15), b(15), c(15), d(15)] = lookup_idw_200n(15);
[w60, x60, y60, z60] = lookup_idw_60n(7);
[w80, x80, y80, z80] = lookup_idw_80n(7);
[a100, b100, c100, d100] = lookup_idw_100n(11);
[a150, b150, c150, d150] = lookup_idw_150n(11);
[w300, x300, y300, z300] = lookup_idw_300n(7);
[w400, x400, y400, z400] = lookup_idw_400n(7);
[a500, b500, c500, d500] = lookup_idw_500n(11);
[w600, x600, y600, z600] = lookup_idw_600n(7);
[w700, x700, y700, z700] = lookup_idw_700n(7);
[w800, x800, y800, z800] = lookup_idw_800n(7);
[w1u, x1u, y1u, z1u] = lookup_idw_1u(7);



ratio = 4;

% gm/Id   
gmid_1 = 11; % 8-12 recommended
gmid_8 = 11; % 8-12 recommened
gmid_nbias = 7; % M3, M6, M7, M9, M13, Mbn1, Mbn2  M14, M17 (M2, M12) % 6-8 recommended
gmid_cmfb = 7; % M4, M10   % 6-8 recommended
gmid_pbias = 7; % M5, Mbp1, M15, maybe M16
gmid_11 = 7;

gm1 = 7e-3;
gm8 = 25e-3;
gm5 = gm1*gmid_pbias/gmid_1;  % gm1*gmid_pbias/gmid_1 means Id1=Id5


Id1 = gm1/gmid_1;
Id8 = gm8/gmid_8;
Id3 = 2*Id1;    gm3 = Id3*gmid_nbias;
Id5 = gm5/gmid_pbias;  
Id4 = Id1+Id5;    gm4 = Id4*gmid_cmfb;
Id6 = Id5; % gm6/gmid_nbias;   
gm6 = Id6*gmid_nbias;  % could be not true so need to choose
Id7=Id6;
Id9 = Id8;
gm7 = Id7*gmid_nbias;
gm9 = Id9*gmid_nbias;
Id13 = Id3;
Id10 = Id13/2;
Id11 = Id10;
gm10 = Id10*gmid_cmfb;
gm11 = Id11*gmid_11;
gm13 = Id13*gmid_nbias;
Idbn1 = Id3/ratio;
Idbn2 = Id3/ratio;
Idbp1 = Id5;
gmbn1 = Idbn1*gmid_nbias;
gmbn2 = Idbn2*gmid_nbias;
gmbp1 = Idbp1*gmid_pbias;
Id14 = Idbp1;
Id15 = Id14;
Id17 = Idbp1;
Id16 = Id15+Idbp1;
gm14 = Id14*gmid_nbias;
gm15 = Id15*gmid_pbias;
gm16 = Id16*gmid_pbias;
gm17 = Id17*gmid_nbias;


W1 = Id1/b100; % b(gmid_1);
W8 = Id8/d150; % d(gmid_8); % d100; % d500;
W3 = Id3/x600; % b(gmid_nbias)
W4 = Id4/d(gmid_cmfb);  % z1; 
W5 = Id5/d(gmid_pbias);
W6 = Id6/b(gmid_nbias); % x400;   % x600; 
W7 = Id7/b(gmid_nbias);
W9 = Id9/x1u;  % b(gmid_nbias);
W10 = Id10/d(gmid_cmfb);
W11 = Id11/b(gmid_11);
W13 = Id13/x600;% b(gmid_nbias);
Wbn1 = Idbn1/x60; % b(gmid_nbias); % x600;% x400; 
Wbn2 = Idbn2/x60; % b(gmid_nbias);
Wbp1 = Idbp1/d(gmid_pbias);
W14 = Id14/x400; % b(gmid_nbias);
W15 = Id15/d(gmid_pbias);
W16 = Id16/d(gmid_pbias);
W17 = Id17/x400; % b(gmid_nbias);

L1 = 100;
L8 = 150;
L3 = 600;
L4 = 200;
L5 = 200;
L6 = 200;
L7 = 200;
L9 = 1000;
L10 = 200;
L11 = 200;
L13 = 600;
Lbn1 = 60; 
Lbn2 = 60;
Lbp1 = 200;
L14 = 400;
L15 = 200;
L16 = 200;
L17 = 400;




% === Build Table ===
Transistor = {'M1';'M3';'M4';'M5';'M6';'M7';'M8';'M9';'M10';'M11';'M13';'Mbn1';'Mbn2'; 'Mbp1'; 'M14'; 'M15'; 'M16'; 'M17'};
gm_ID = [gmid_1; gmid_nbias; gmid_cmfb; gmid_pbias; gmid_nbias; gmid_nbias; ...
         gmid_8; gmid_nbias; gmid_cmfb; gmid_11; gmid_nbias; gmid_nbias; gmid_nbias; gmid_pbias;...
         gmid_nbias; gmid_pbias; gmid_pbias; gmid_nbias];

gm_values = [gm1; gm3; gm4; gm5; gm6; gm7; gm8; gm9; gm10; gm11; gm13; gmbn1; gmbn2; gmbp1;...
    gm14; gm15; gm16; gm17];

Id_values = [Id1; Id3; Id4; Id5; Id6; Id7; Id8; Id9; Id10; Id11; Id13; Idbn1; Idbn2; Idbp1;...
    Id14; Id15; Id16; Id17];

W_values  = [W1; W3; W4; W5; W6; W7; W8; W9; W10; W11; W13; Wbn1; Wbn2; Wbp1; W14; W15; W16; W17];

L_values_nm  = [L1; L3; L4; L5; L6; L7; L8; L9; L10; L11; L13; Lbn1; Lbn2; Lbp1; L14; L15; L16; L17];

T = table(Transistor, gm_ID, gm_values, Id_values, W_values, L_values_nm)

% === Save to CSV (optional) ===
% writetable(T, 'transistor_sizing_table.csv');




