function truth= gen_truth(model)

%variables
truth.K= 100;                   %length of data/number of scans
truth.X= cell(truth.K,1);             %ground truth for states of targets  
truth.N= zeros(truth.K,1);            %ground truth for number of targets
truth.L= cell(truth.K,1);             %ground truth for labels of targets (k,i)
truth.track_list= cell(truth.K,1);    %absolute index target identities (plotting)
truth.total_tracks= 0;          %total number of appearing tracks

%target initial states and birth/death times
nbirths= 1;

xstart(:,1)  = [ 0; 1; 0; 7; 0.03 ];        tbirth(1)  = 10;     tdeath(1)  = 80;

%generate the tracks
for targetnum=1:nbirths
    targetstate = xstart(:,targetnum);
    for k=tbirth(targetnum):min(tdeath(targetnum),truth.K)
        targetstate = gen_newstate_fn(model,targetstate,'noiseless');
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
     end
end
truth.total_tracks= nbirths;

