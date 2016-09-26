
NumberOfSamples = 10;
NumberOfElements = 26640;

cols = [];
rows = [];
values = [];

file = fopen('C:/temp/triplets.txt', 'r');

shape = [3 Inf];
formatSpec = '%d %d %f';

disp('Loading samples position info...');
fromFile = fscanf(file, formatSpec, shape);
    
rows = fromFile(1,:) + 1;
cols = fromFile(2,:) + 1;
values = fromFile(3,:);

rows = [rows, max(rows)];
%cols = [cols, 26460];
cols = [cols, 3087];

values = [values, 0];

disp('Building A, At and AtA...');

A = sparse(rows, cols, values); 
At = transpose(A);
AtA = At * A;

counter = 1;

for channel = 0:2
    b = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    file = fopen(strcat('C:/temp/b',num2str(channel),'.txt'), 'r');

%     tline = fgetl(file);
    disp(strcat('Loading samples irradiance info (channel ',num2str(channel),')...'));
    shape = [1 5];
    formatSpec = '%f';
    b = fscanf(file, formatSpec);
%     while ischar(tline)
%         %disp(tline)
%         C = strsplit(tline);
% 
%         value = str2double(C{1});
%         b = [b, value];
% 
% 
%         tline = fgetl(file);
%         counter = counter + 1;
%     end

    fclose(file);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    guess = [];
    file = fopen(strcat('C:/temp/guess',num2str(channel),'.txt'), 'r');

    tline = fgetl(file);
    disp(strcat('Loading 1st guess (channel ',num2str(channel),')...'));

    while ischar(tline)
        %disp(tline)
        C = strsplit(tline);

        value = str2double(C{1});
        guess = [guess, value];


        tline = fgetl(file);
    end

    fclose(file);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Atb = At * b;

    disp(strcat('Attempting to linear solve for channel ', num2str(channel),'...'));

    tol = 1e-6;
    maxit = 1e7;
    %L = ichol(A);
    [x,flag,err,iter,res] = bicgstab(AtA,Atb,tol,maxit, [], [], transpose(guess));
    %[x,flag,err,iter,res] = bicgstab(AtA,Atb,tol,maxit);

    dlmwrite(strcat('C:/temp/v',num2str(channel),'.txt'), x, 'delimiter', '\n');
    %save(strcat('C:/temp/v',num2str(channel),'.txt'), 'x', '-ascii', '-double');
end