function r = call_on(func, array)
  t = num2cell(array);
  r = func(t{:});
end