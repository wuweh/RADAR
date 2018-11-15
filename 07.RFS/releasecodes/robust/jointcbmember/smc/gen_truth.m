function truth= gen_truth(model)

%variables
truth.K= 100;                   %length of data/number of scans
truth.X= cell(truth.K,1);             %ground truth for states of targets  
truth.N= zeros(truth.K,1);            %ground truth for number of targets
truth.L= cell(truth.K,1);             %ground truth for labels of targets (k,i)
truth.track_list= cell(truth.K,1);    %absolute index target identities (plotting)
truth.total_tracks= 0;          %total number of appearing tracks

%target initial states and birth/death times
nbirths= 10;
sigma_jiggle = 0;
u_jiggle= 100; v_jiggle= 2;
wturn = 2*pi/180;

xstart(:,1)  = [0.98; 1000+3.8676; -10; 1500-11.7457; -10; wturn/8 ];        tbirth(1)  = 1;     tdeath(1)  = truth.K+1;
xstart(:,2)  = [0.98; -250-5.8857;  20; 1000+11.4102; 3; -wturn/3 ];         tbirth(2)  = 10;    tdeath(2)  = truth.K+1;
xstart(:,3)  = [0.98; -1500-7.3806; 11; 250+6.7993; 10; -wturn/2 ];          tbirth(3)  = 10;    tdeath(3)  = truth.K+1;
xstart(:,4)  = [0.98; -1500; 43; 250; 0; 0 ];                                tbirth(4)  = 10;    tdeath(4)  = 66;
xstart(:,5)  = [0.98; 250-3.8676; 11; 750-11.0747; 5; wturn/4 ];             tbirth(5)  = 20;    tdeath(5)  = 80;
xstart(:,6)  = [0.98; -250+7.3806; -12; 1000-6.7993; -12; wturn/2 ];         tbirth(6)  = 40;    tdeath(6)  = truth.K+1;
xstart(:,7)  = [0.98; 1000; 0; 1500; -10; wturn/4 ];                         tbirth(7)  = 40;    tdeath(7)  = truth.K+1;
xstart(:,8)  = [0.98; 250; -50; 750; 0; -wturn/4 ];                          tbirth(8)  = 40;    tdeath(8)  = 80;
xstart(:,9)  = [0.98; 1000; -50; 1500; 0; -wturn/4 ];                        tbirth(9)  = 60;     tdeath(9)  = truth.K+1;
xstart(:,10)  = [0.98; 250; -40; 750; 25; wturn/4 ];                         tbirth(10)  = 60;    tdeath(10)  = truth.K+1;

%generate the tracks
for targetnum=1:nbirths
    targetstate = xstart(:,targetnum);
    for k=tbirth(targetnum):min(tdeath(targetnum),truth.K)
        exc= betarnd(u_jiggle,v_jiggle,1)-targetstate(1); exc=0;
        excitation= [zeros(2,1); sigma_jiggle*randn(1)];
        targetstate = gen_newstate_fn_tg(model,targetstate,exc,excitation);
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
     end
end
truth.total_tracks= nbirths;


