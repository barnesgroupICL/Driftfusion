function moviewrite(Framefile, name)
% Write a frame file to an avi movie
% name is a string with the desired filename- do NOT include .avi extension

% name is a string for the final video
name = [name, '.avi'];

% Write to file
myVideo = VideoWriter(name);
myVideo.FrameRate = 6;
open(myVideo);
writeVideo(myVideo, Framefile);
close(myVideo);

end