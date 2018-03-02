function img2mp4( img, framerate,videoname)
% convert images to mp4 video
% input :
%     img: Nx*Ny*1*Nt or Nx*Ny*Nt for grayscale or Nx*Ny*3*Nt for RGB
%     framerate: frames per second

%    Ziyue Wu 04/03/2014

if length(size(img)) == 3
    img = permute(img, [1 2 4 3]);
end
writerObj = VideoWriter(videoname,'MPEG-4');
writerObj.FrameRate = framerate;
open(writerObj);
writeVideo(writerObj,img);
close(writerObj);
end