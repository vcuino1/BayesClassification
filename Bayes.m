%% ���ݷ��࣬��Ϊѵ�����Ͳ��Լ�
clear;
clc;
temp = 0;
traindata = load('traindata.mat');
traindata = traindata.traindata;
testdata = load('testdata.mat');
testdata = testdata.testdata;
%% �����i������ȡjֵ�ĸ��ʣ�CN1�Ƿ���Ϊ��һ��ǰ����ȡֵ��������������
for i = 1:4
    A = traindata(1:45,i);
    B = unique(A);
    C= zeros(size(B));
    for j = 1:length(B)
        C(j) = length(find(A==B(j)));
        CN1(j,2*i-1) = B(j);
        CN1(j,2*i) = C(j)/45;
    end
end
for i = 1:4
    A = traindata(46:90,i);
    B = unique(A);
    C= zeros(size(B));
    for j = 1:length(B)
        C(j) = length(find(A==B(j)));
        CN2(j,2*i-1) = B(j);
        CN2(j,2*i) = C(j)/45;
    end
end
for i = 1:4
    A = traindata(91:135,i);
    B = unique(A);
    C= zeros(size(B));
    for j = 1:length(B)
        C(j) = length(find(A==B(j)));
        CN3(j,2*i-1) = B(j);
        CN3(j,2*i) = C(j)/45;
    end
end
%% ����
for i = 1:15
    Statistic1(i,:) = 1;
    for j = 1:4
        if isempty(find(CN1(:,2*j-1) == testdata(i,j)))
            [~,a] = min(abs(CN1(:,2*j-1)-testdata(i,j)));%������Լ���ĳ����ֵ��ѵ�����в����ڣ��Ǿ���ѵ������Ѱ�����ֵ��ӽ�����ֵ
            Position1(i,j) = a;%min(find(CN1(:,2*j-1) == testdata(i,j)));
        else
            Position1(i,j) = find(CN1(:,2*j-1) == testdata(i,j));
        end
        Statistic1(i,:) = CN1(Position1(i,j),2*j)*Statistic1(i,:);
    end
end
for i = 1:15
    Statistic2(i,:) = 1;
    for j = 1:4
        if isempty(find(CN2(:,2*j-1) == testdata(i,j)))
            [~,a] = min(abs(CN2(:,2*j-1)-testdata(i,j)));%������Լ���ĳ����ֵ��ѵ�����в����ڣ��Ǿ���ѵ������Ѱ�����ֵ��ӽ�����ֵ
            Position2(i,j) = a;%min(find(CN2(:,2*j-1) == testdata(i,j)));
        else
            Position2(i,j) = find(CN2(:,2*j-1) == testdata(i,j));
        end
        Statistic2(i,:) = CN2(Position2(i,j),2*j)*Statistic2(i,:);
    end
end
for i = 1:15
    Statistic3(i,:) = 1;
    for j = 1:4
        if isempty(find(CN3(:,2*j-1) == testdata(i,j)))
            [~,a] = min(abs(CN3(:,2*j-1)-testdata(i,j)));%������Լ���ĳ����ֵ��ѵ�����в����ڣ��Ǿ���ѵ������Ѱ�����ֵ��ӽ�����ֵ
            Position3(i,j) = a;%min(find(CN3(:,2*j-1) == testdata(i,j)));
        else
            a=find(CN3(:,2*j-1) == testdata(i,j));
            Position3(i,j) = min(find(CN3(:,2*j-1) == testdata(i,j)));
        end
        Statistic3(i,:) = CN3(Position3(i,j),2*j)*Statistic3(i,:);
    end
end
% �ж�
for i = 1:15
    if Statistic1(i)>Statistic2(i) && Statistic1(i)>Statistic3(i)
        label(i,:) = 1;
        if label(i,:) == testdata(i,5)
            temp = temp + 1;
        end
    elseif Statistic2(i)>Statistic1(i) && Statistic2(i)>Statistic3(i)
        label(i,:) = 2;
        if label(i,:) == testdata(i,5)
            temp = temp + 1;
        end
    else Statistic3(i)>Statistic1(i) && Statistic3(i)>Statistic2(i)
        label(i,:) = 3;
        if label(i,:) == testdata(i,5)
            temp = temp + 1;
        end
    end
end
rate = temp/15;