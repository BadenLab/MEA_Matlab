%Calculate the beginning of a stimulus relative to the beginning of the
%trace
function  Out_times = Find_Epoch_trace (epochs, tbegin, tend, binsize)
a = 1;

for i = epochs
    Trace_begin(a) = ceil(tbegin(i)-tbegin(1));
    Trace_end(a) = ceil(tend(i)-tbegin(1));
    a = a+1;
    
end

%Calculate the location in the binned trace
Trace_begin_b = Trace_begin/binsize;
Trace_end_b = Trace_end/binsize;

Out_times.Trace_begin = Trace_begin;
Out_times.Trace_end = Trace_end;
Out_times.Trace_begin_b = Trace_begin_b;
Out_times.Trace_end_b = Trace_end_b;


      
    
end
