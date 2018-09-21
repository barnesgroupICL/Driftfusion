function writemovie(Framefile, name)

% name is a string for the final video
name = [name, '.avi'];

% Write to file
myVideo = VideoWriter(name);
myVideo.FrameRate = 10;
open(myVideo);
writeVideo(myVideo, Framefile);
close(myVideo);

end