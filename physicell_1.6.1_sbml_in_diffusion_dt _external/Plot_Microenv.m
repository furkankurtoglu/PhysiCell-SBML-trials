close all
clear
clc

cd output

s=what;
MatFiles = s.mat;
OutMatFiles = MatFiles(contains(MatFiles,'micro'));
OutMatFiles(1) = [];
OutMatFiles(1) = [];
for i = 1:length(OutMatFiles)
    OutMatFiles{i}=OutMatFiles{i}(1:14);
end


%%
for i = 1:length(OutMatFiles)

    xmlname=strcat(OutMatFiles{i},'.xml');
    MCDS = read_MultiCellDS_xml( xmlname);
    
    %% Oxygen
    k = find( MCDS.mesh.Z_coordinates == 0 );
    figure(1)
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), MCDS.continuum_variables(1).data(:,:,k) ,20 ) ;
    axis image
    colorbar
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) );
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) );

    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(1).name , ...
        MCDS.continuum_variables(1).units , ...
        MCDS.metadata.current_time , ...
        MCDS.metadata.time_units, ...
        MCDS.mesh.Z_coordinates(k), ...
        MCDS.metadata.spatial_units ) );

    
    %% Glucose
    k = find( MCDS.mesh.Z_coordinates == 0 );
    figure(2)
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), MCDS.continuum_variables(2).data(:,:,k) ,20 ) ;
    axis image
    colorbar
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) );
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) );

    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(2).name , ...
        MCDS.continuum_variables(2).units , ...
        MCDS.metadata.current_time , ...
        MCDS.metadata.time_units, ...
        MCDS.mesh.Z_coordinates(k), ...
        MCDS.metadata.spatial_units ) );
    

    %% Lactate
    k = find( MCDS.mesh.Z_coordinates == 0 );
    figure(3)
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), MCDS.continuum_variables(3).data(:,:,k) ,20 ) ;
    axis image
    colorbar
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) );
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) );

    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(3).name , ...
        MCDS.continuum_variables(3).units , ...
        MCDS.metadata.current_time , ...
        MCDS.metadata.time_units, ...
        MCDS.mesh.Z_coordinates(k), ...
        MCDS.metadata.spatial_units ) );
    
    
    
end 

cd ..