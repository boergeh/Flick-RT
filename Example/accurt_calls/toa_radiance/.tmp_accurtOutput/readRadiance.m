function data = readRadiance(fileName)

  fid = fopen(fileName);

  % Number of repeated runs (see the REPEATED_RUN_SIZE tag in the main
  % configuration file).
  nRuns = fscanf(fid,'%i',1);

  for i = 1:nRuns
    data(i).nStreams       = fscanf(fid,'%i',1);
    data(i).nDepths        = fscanf(fid,'%i',1); 
    data(i).nWavelengths   = fscanf(fid,'%i',1); 
    data(i).nPolarAngles   = fscanf(fid,'%i',1); 
    data(i).nAzimuthAngles = fscanf(fid,'%i',1); 

    
    for j = 1:data(i).nDepths
      data(i).depths(j) = fscanf(fid,'%g',1); % [m]
    end
    for j = 1:data(i).nWavelengths
      data(i).wavelengths(j) = fscanf(fid,'%g',1); % [nm]
    end
    for j = 1:data(i).nPolarAngles
      data(i).polarAngles(j) = fscanf(fid,'%g',1); % [degrees]
    end
    for j = 1:data(i).nAzimuthAngles
      data(i).azimuthAngles(j) = fscanf(fid,'%g',1); % [degrees]
    end
    
    for j = 1:data(i).nDepths
      for k = 1:data(i).nWavelengths
	for l = 1:data(i).nPolarAngles
	  for m = 1:data(i).nAzimuthAngles
	    data(i).radiance(j,k,l,m) = fscanf(fid,'%g',1); % [W/m2/nm/sr]
	  end
	end
      end
    end
  end

  fclose(fid);


  
  