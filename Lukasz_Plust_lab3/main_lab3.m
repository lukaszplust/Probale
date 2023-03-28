clc
clear all
close all

% odpowiednie fragmenty kodu mozna wykonac poprzez zaznaczenie i wcisniecie F9 w Matlabie
% komentowanie/odkomentowywanie: ctrl+r / ctrl+t
%na macu fn+shift+f7
% Zadanie A
%------------------
N = 10; %liczna stron
density = 3; % parametr decydujacy o gestosci polaczen miedzy stronami
[Edges] = generate_network(N, density);
%-----------------

% Zadanie B
%------------------
% generacja macierzy I, A, B i wektora b
% macierze A, B i I musza byc przechowywane w formacie sparse (rzadkim)

d = 0.85;
B = sparse(Edges(2,:),Edges(1,:),1,N,N);
I = speye(N);
L = sparse(sum(B))';% zmienilem L z poprzedniego zadania w celu ulatwienia operacji w petli
b = zeros(N,1);
b(:,1) = (1-d)/N;
A=spdiags(1./L,0,N,N);
M = sparse(I - d*B*A);

%issparse(M);

%-----------------
%Zadanie C
%------------------
r = M \ b;

%------------------

% Zadanie D
%------------------
close all

N = [500, 1000, 3000, 6000, 12000];
density = 10;

for i = 1:5

    [Edges] = generate_network(N(i),density);
    B = sparse(Edges(2,:),Edges(1,:),1,N(i),N(i));
    I = speye(N(i));
    L = sparse(sum(B))';
    b = zeros(N(i),1);
    b(:,1) = (1-d)/N(i);
    A=spdiags(1./L,0,N(i),N(i));
    M = sparse(I - d*B*A);

    tic
    % obliczenia start
    r = M \ b; %pomiar czasu obejmuje tylko rozwiazanie rownania
    % obliczenia stop
    czas_Gauss(i) = toc;
end

figure("Name", "Time Gauss");
plot(N, czas_Gauss)
title("Time solution of Gauss method");
xlabel("Size of martix N");
ylabel("Time in (s)");
saveas(gcf,"zadD.png");
%------------------



% Zadanie E
%------------------
clc
close all

% sprawdz przykladowe dzialanie funkcji tril, triu, diag:

%Z = rand(4,4);%losowanie wartosci
%tril(Z,-1);%dolna trojkatna
%triu(Z,1); %gorna trojkatna
%diag(diag(Z));%wartosci na diagonali

%N = [500, 1000, 3000, 6000, 12000]; %dodac jakbym odpalal bez wczesniej zadeklarowanego
%density = 10; %dodac jakbym odpalal bez wczesniej zadeklarowanego
%d = 0.85; %dodac jakbym odpalal bez wczesniej zadeklarowanego

condition = 10^(-14);

for i = 1:5
    
    [Edges] = generate_network(N(i),density);
    B = sparse(Edges(2,:),Edges(1,:),1,N(i),N(i));
    I = speye(N(i));
    L = sparse(sum(B))';
    b = zeros(N(i),1);
    b(:,1) = (1-d)/N(i);
    A=spdiags(1./L,0,N(i),N(i));
    M = sparse(I - d*B*A);

    tic
    % obliczenia start
    
    r = ones(N(i),1);% inicjalizacja wektora r
    D = diag(diag(M)); %macierz złożona tylko z elementów diagonali
    L = tril(M,-1); %macierz trojkatna dolna (bez diagonali jak z pdf)
    U = triu(M,1); % macierz trojkatna gorna (bez diagonali jak z pdf)
    
    iterations(i)= 0;%licznik iteracji

    while(true) %dopoki norma bledu rezydualnego jest mniejsza od 10^-14
        iterations(i) = iterations(i) +1; %zwiekszam licznik operacji
        r = (-D \(L + U))*r + (D \ b); %wyznaczam ze wzoru z pdf
        res = M * r -b; % wzor z pdf
        %przypisuje do residual_norm_Jacobi wartosci norm(res) dla kazdej
        %iteracji
        residual_norm_Jacobi(iterations(i)) = norm(res);
        if norm(res) < condition
            break;
        end
    end 
    % obliczenia stop
    czas_Jacobi(i) = toc;
end

%wykres czasu
figure("Name", "Time Jacobi");
plot(N, czas_Jacobi)
title("Time solution of Jacobi method");
xlabel("Size of martix N");
ylabel("Time in [s]");
saveas(gcf,"zadE_czas.png");

%wykres iteracji
figure("Name", "Iterations of Jacobi method");
plot(N, iterations);%216
title("Number of iterations of Jacobi method");
xlabel("Size of martix N");
ylabel("Number of iterations");
saveas(gcf,"zadE_iteracje.png");

%norma
figure("Name", "Residual norm of Jacobi method");
% iterations 207   209   212   214   216
% wyswietlam residual_norm_Jacobi od 1:216
semilogy(residual_norm_Jacobi(1:iterations(i))); 
title("Residual norm of Jacobi method for N=1000");
xlabel("Number of iterations");
ylabel("Residual norm")
saveas(gcf,"zadE_norma.png");


%------------------


% Zadanie F
%------------------
close all
clc
condition = 10^(-14); %-> dodac jesli nie było dodane

for i = 1:5
    
    [Edges] = generate_network(N(i),density);
    B = sparse(Edges(2,:),Edges(1,:),1,N(i),N(i));
    I = speye(N(i));
    L = sparse(sum(B))';
    b = zeros(N(i),1);
    b(:,1) = (1-d)/N(i);
    A=spdiags(1./L,0,N(i),N(i));
    M = sparse(I - d*B*A);

    tic
    % obliczenia start
    
    r = ones(N(i),1);% inicjalizacja wektora r
    D = diag(diag(M)); %macierz złożona tylko z elementów diagonali
    L = tril(M,-1); %macierz trojkatna dolna (bez diagonali jak z pdf)
    U = triu(M,1); % macierz trojkatna gorna (bez diagonali jak z pdf)
    
    iterations(i)= 0;%licznik iteracji

    while(true) %dopoki norma bledu rezydualnego jest mniejsza od 10^-14
        iterations(i) = iterations(i) +1; %zwiekszam licznik operacji
        r = (-(D + L))\(U*r) + ((D + L) \ b); % 7 pdpkt. wyznaczam ze wzoru z pdf
        res = M * r -b; % wzor z pdf
        %przypisuje do residual_norm_Gauss_Seidel wartosci norm(res) dla kazdej
        %iteracji
        residual_norm_Gauss_Seidel(iterations(i))= norm(res);
        if norm(res) < condition
            break;
        end
    end 
    % obliczenia stop
    czas_Gauss_Seidel(i) = toc;
end

%wykres czasu
figure("Name", "Time Gauss Seidel");
plot(N, czas_Gauss_Seidel)
title("Time solution of Gauss Seidel method");
xlabel("Size of martix N");
ylabel("Time in [s]");
saveas(gcf,"zadF_czas.png");

%wykres iteracji
figure("Name", "Iterations of Gauss Seidel method");
plot(N, iterations);%112
title("Number of iterations of Gauss Seidel method");
xlabel("Size of martix N");
ylabel("Number of iterations");
saveas(gcf,"zadF_iteracje.png");

%norma
figure("Name", "Residual norm of Gauss Seidel method ");
% iterations = 107   109   110   111   112
%wyswietlam residual_norm_Gauss_Seidel od 1:112
semilogy(residual_norm_Gauss_Seidel(1:iterations(i)));
title("Residual norm of Gauss Seidel method for N=1000");
xlabel("Number of iterations");
ylabel("Residual norm")
saveas(gcf,"zadF_norma.png");

%------------------

% Zadanie G
%------------------
load("Dane_Filtr_Dielektryczny_lab3_MN.mat");

%metoda Gaussa
r = M \ b;
save result_Gauss r


%------------------
%metoda Jacobiego

condition = 10^(-14);

D = diag(diag(M)); 
U = triu(M, 1); 
L = tril(M, -1);

%to samo co wyzej tylko bez iteracji :)
while(true)
    r = (-D \(L + U))*r + (D \ b); %wyznaczam ze wzoru z pdf
    res = M * r -b; % wzor z pdf
    if norm(res) < condition
        break;
    end
end

save result_Jacobi r

% Norma idzie do nieskonczoności az zamieni sie w Nan.

%------------------
%metoda Gaussa-Seidla

%to dodac jesli odpalam osobno
%------------------
%condition = 10^(-14);

%D = diag(diag(M)); 
%U = triu(M, 1); 
%L = tril(M, -1);
%------------------
while(true) %dopoki norma bledu rezydualnego jest mniejsza od 10^-14
    r = (-(D + L))\(U*r) + ((D + L) \ b); % 7 pdpkt. wyznaczam ze wzoru z pdf
    res = M * r -b; % wzor z pdf
    if norm(res) < condition
        break;
    end
end

save result_Gauss_Seidel r

% Norma idzie do nieskonczoności az zamieni sie w Nan.