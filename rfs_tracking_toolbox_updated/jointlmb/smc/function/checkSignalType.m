function [ signalType ] = checkSignalType( Tp,f0,type )
Tps = [100,200,400,800,1600];
f0s = [2.73 3 3.28];
types = [string('FM3'),string('CW3'),string('FM'),string('CW')];
Tp = find(Tps==Tp);f0 = find(f0s==f0);type = find(types==type);
if isempty(Tp)||isempty(f0)||isempty(type)
    signalType = -1;
else
    switch Tp*100+f0*10+type
        case 311
            signalType = 1;
        case 331
            signalType = 2;
        case 322
            signalType = 3;
        case 321
            signalType = 4;
        case 323
            signalType = 5;
        otherwise
            disp(Tp*100+f0*10+type)
            signalType = -1;
    end  
end

