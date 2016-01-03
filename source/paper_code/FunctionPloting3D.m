%Omid55
c = 1;
x1 = -5:0.1:5;
x2 = -5:0.1:5;
d1 = zeros(length(x1)*length(x2),1);
d2 = zeros(length(x1)*length(x2),1);
d3 = zeros(length(x1)*length(x2),1);
mat = zeros(length(x1),length(x2));
for i=1:length(x1)
    for j=1:length(x2)
        d1(c) = x1(i);
        d2(c) = x2(j);
        r = MyFunctionForPloting3D(x1(i)) + MyFunctionForPloting3D(x2(j));
        mat(i,j) = r;
        d3(c) = r;
        c = c + 1;
    end
end
figure, plot3(d1,d2,d3);
figure, surf(x1,x2,mat);