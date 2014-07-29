 normalized_throughput = [];
for (delta = 1:1) %1:5
for (loc_no= 1:1) %1:100
    delta
    loc_no
  normalized_throughput(delta,loc_no)  =   maxthrough_directional_deter(delta,loc_no);
end
end