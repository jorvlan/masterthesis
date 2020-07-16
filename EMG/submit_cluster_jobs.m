
addpath '/home/common/matlab/fieldtrip/qsub'
cd /project/3024005.01/Analysis/Tremor/log


subjects = {'002','003','004','005','006','007','008','010','011','012','013','014','015','016','017','018'};


jobs = {};
req_mem   = 10^10; % 10 GB
req_etime = 7200; % 2 hours - hopefully scheduled slightly faster
mscript = '/project/3024005.01/Analysis/Tremor/tremor_1st_level_TOS';


for i = 1:length(subjects)
    jobs{i} = qsubfeval(mscript,  subjects{i}, 'memreq', req_mem, 'timreq', req_etime);
end

save '/project/3024005.01/Analysis/Tremor/log/jobs.mat' jobs
