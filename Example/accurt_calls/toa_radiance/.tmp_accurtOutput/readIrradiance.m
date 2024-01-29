function data = readIrradiance(fileName)

  fid = fopen(fileName);

  % Number of repeated runs (see the REPEATED_RUN_SIZE tag in the main
  % configuration file).
  nRuns = fscanf(fid,'%i',1);
    
  for i = 1:nRuns
    data(i).nStreams     = fscanf(fid,'%i',1);
    data(i).nDepths      = fscanf(fid,'%i',1); 
    data(i).nWavelengths = fscanf(fid,'%i',1); 
    
    for j = 1:data(i).nDepths
      data(i).depth(j) = fscanf(fid,'%g',1); % [m]
    end
    for j = 1:data(i).nWavelengths
      data(i).wavelength(j) = fscanf(fid,'%g',1); % [nm] 
    end
    
    for j = 1:data(i).nDepths
      for k = 1:data(i).nWavelengths
	data(i).irradiance(j,k) = fscanf(fid,'%g',1); % [W/m2/nm]
      end
    end
  end

  fclose(fid);


  
  