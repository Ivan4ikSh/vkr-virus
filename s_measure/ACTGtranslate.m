%Made with MATLAB2022a
%% parameters
% b1 - number of the first sequence to be read
% b2 - number of the last sequnce to be read
% tr1 - site gaps threshold
% tr2 - sequence gaps threshold
% fcut - global monomorphous threshhold
% fcutloc - local monomorphous threshhold
b1 = 1; 
b2 = 100;
tr1 = 0.05; 
tr2 = 0.1;
fcut = 0.1;
fcutloc = 0.02;
lenmode = true; %mode of sequence different length regulation
gen = {}; %container for raw info
data = {}; %container for result info
l = []; %container for sequence length
%take data from files
 
for k = 1:3
    am = (b2-b1+1); %amount of the sequences
    filepath=strcat('data\period',num2str(k));
    oldpath=cd(filepath);
filelist{k}= dir ('**/*.fas');
gen{k} = fastaread(filelist{k}.name,'BlockRead',[b1 b2],'IgnoreGaps',false);
 l = [l length(gen{k}(1).Sequence)];
%l - length of the sequence
cd(oldpath);
end
disp(tr1*am*3) 
% sequence length regulation
% true - cutoff
% false - extend
if(lenmode)
    l = min(l);
else
for k = 1:3
    for l1 = 1:(max(l) - l(k))
        for j = 1:am
            gen{k}(j).Sequence = convertStringsToChars(gen{k}(j).Sequence + "-");
        end
    end

end
l = max(l)
end



    %data cleaning
    banlistSeq = [];
banlistSit = [];

for i0 = 1:l
    del = 0;
for k = 1:3
for j0 = 1:am
    seq = gen{k}(j0).Sequence;

if seq(i0) ~='A' &&  seq(i0) ~='C' && seq(i0) ~='T' && seq(i0) ~='G'

    del = del+1;
end
end
end

if del ~=0
 if del/(am*3) >=tr1

banlistSit = [banlistSit i0];

end

end

end


if ~isempty(banlistSit)
    for k = 1:3
    for sch3 = 1:am
  gen{k}(sch3).Sequence(banlistSit) = 'E';  
    end
    end
end


for k = 1:3
   am = (b2-b1+1);
for j1 = 1:am
    del = 0;
    for i1 = 1:l
        seq = gen{k}(j1).Sequence;
if seq(i1) ~='A' &&  seq(i1) ~='C' && seq(i1) ~='T' && seq(i1) ~='G' && seq(i1) ~='E'
    del = del+1;
end
    end

if del ~=0

 if del/l >=tr2
banlistSeq = [banlistSeq j1];
 end
end

end


if ~isempty(banlistSeq)

gen{k}(banlistSeq) = [];
am = am - length(banlistSeq);
banlistSeq = [];
end
data{k} = zeros(am,l);
end
disp(tr2*l)


%finding the common consensus 
cons = '';
for i = 1:l
    qA=0;
    qC=0;
    qT=0;
    qG=0;
    qE=0;
    for k = 1:3
    for j = 1:length(gen{k})
        seq = gen{k}(j).Sequence;
        switch seq(i)
            case 'A'
qA = qA+1;
            case 'C'
qC = qC+1;
            case 'T'
qT = qT+1;
            case 'G'
qG = qG+1;    
            case 'E'
qE = qE+1;
        end
    end

    end

    [~,I] = max([qA qC qT qG qE]);
switch I
    case 1
cons = [cons 'A'];
    case 2
cons = [cons 'C'];        
    case 3
cons = [cons 'T'];        
    case 4
cons = [cons 'G'];
    case 5
cons = [cons 'E'];
end
end



%binarization
for i2 = 1:l
for k = 1:3
    for j2 = 1:length(gen{k})
        seq = gen{k}(j2).Sequence;
if seq(i2)~='E'
        if seq(i2) == cons(i2)
            data{k}(j2,i2)=0;
        else
            data{k}(j2,i2)=1;
        end
else
    data{k}(j2,i2)=2;
end
    end
end
end



%deliting the monomorph and marked "bad" sites
rightsites=[];
disp(fcut*(length(data{1}(:,1))+length(data{2}(:,1))+length(data{3}(:,1))))

for msch = 1:l
if data{1}(1,msch)~=2
    %if ~(mean(data{k}(:,msch))>=1-fcut || mean(data{k}(:,msch))<=fcut)
    av = mean([mean(data{1}(:,msch)), mean(data{2}(:,msch)), mean(data{3}(:,msch))]);
    if ~(av<=fcut || av>=1-fcut)        
rightsites = [rightsites msch];
    end
end
end
disp([fcutloc*length(data{1}(:,1)), (1-fcutloc)*length(data{1}(:,1))])
disp([fcutloc*length(data{2}(:,1)), (1-fcutloc)*length(data{2}(:,1))])
disp([fcutloc*length(data{3}(:,1)), (1-fcutloc)*length(data{3}(:,1))])
for k = 1:3
    rightsites2{k}=[];
    for msch = rightsites
            av = mean(data{k}(:,msch));
    if ~(av<=fcutloc || av>=1-fcutloc)        
rightsites2{k} = [rightsites2{k} msch];
    end
    end
end

    %findind the intersection between 'right' sites of the different time
    %periods
    C1=intersect(rightsites2{1},rightsites2{2});
C2=intersect(C1,rightsites2{3});
%removing 'bad' sites and remebering the numbers of the 'right' one
for kk=1:k
  %  [~,l]=size(data{kk});
    C3 = 1:l;
    C3(C2) = [];
    data{kk}(:,C3) = [];
end




data{k+1} =C2;
%saving the data
save(strcat('BinData','.mat'),'data');
disp(data)