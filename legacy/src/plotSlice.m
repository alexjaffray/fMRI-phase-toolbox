function plotSlice(phas,uphas,fl,harmfields,selectedSlice)

    % Plot 3 Images Side by Side
    figure();
    subplot(2,2,1);
    imagesc(phas(:,:,selectedSlice,5));
    xlabel('x-position [voxels]');
    ylabel('y-position [voxels]');
    subplot(2,2,2);
    imagesc(uphas(:,:,selectedSlice,5));colormap(gray)
    xlabel('x-position [voxels]');
    ylabel('y-position [voxels]');
    subplot(2,2,3);
    imagesc(fl(:,:,selectedSlice,5));colormap(gray)
    xlabel('x-position [voxels]');
    ylabel('y-position [voxels]');
    subplot(2,2,4);
    imagesc(harmfields(:,:,selectedSlice,5));colormap(gray)
    xlabel('x-position [voxels]');
    ylabel('y-position [voxels]');

end

