function plot_compress(x,y,col)

for i=1:size(y,2)
    I = find(diff(y(:,i))); I = sort([1; I(:);I(:)+1; size(y,1)]);
    plot(x(I),y(I,i),'color', col)    
end

end