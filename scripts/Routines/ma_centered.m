function y = ma_centered(x, h)
%MA Centered moving average

n = length(x);
y = NaN(size(x));
for i = 1:n
    try
       y(i) = mean(x(i-h:i+h), 'omitnan');
    end
end

end

