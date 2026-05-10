function [ output ] = maxshow( input,X,Y )
fun = @(block_struct) max(block_struct.data(:));
output = blockproc(input,[X,Y],fun);
end