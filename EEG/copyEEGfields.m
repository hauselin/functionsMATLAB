function [copyto] = copyEEGfields(copyfrom,copyto)
   
copyto.subject = copyfrom.subject;
copyto.nbchan = copyfrom.nbchan;
copyto.times = copyfrom.times;
copyto.xmin = copyfrom.xmin;
copyto.xmax = copyfrom.xmax;
copyto.srate = copyfrom.srate;
copyto.pnts = copyfrom.pnts;
copyto.pnts = copyfrom.pnts;
copyto.chanlocs = copyfrom.chanlocs;
copyto.chaninfo = copyfrom.chaninfo;
copyto.ref = copyfrom.ref;
copyto.trials = copyfrom.trials;
copyto.designmat = copyfrom.designmat;

end

