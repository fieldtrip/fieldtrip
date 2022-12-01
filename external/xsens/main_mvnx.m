% load data
tree = load_mvnx(filename);
% read some basic data from the file
mvnxVersion = tree;
fileComments = tree.subject.comment;
%read some basic properties of the subject;
frameRate = tree.subject.frameRate;
suitLabel = tree.subject.label;
originalFilename = tree.subject.originalFilename;
recDate = tree.subject.recDate;
segmentCount = tree.subject.segmentCount;
%retrieve sensor labels
%creates a struct with sensor data
if isfield(tree.subject,'sensors') && isstruct(tree.subject.sensors)
    sensorData = tree.subject.sensors.sensor;
end
%retrieve segment labels
%creates a struct with segment definitions
if isfield(tree.subject,'segments') && isstruct(tree.subject.segments)
    segmentData = tree.subject.segments.segment;
end
%retrieve the data frames from the subject
nSamples = length(tree.subject.frames.frame);
%pre allocate
%pre allocate some memory for the position of Segment1
p_Segment1 = zeros(nSamples,3);
%read the data from the structure e.g. segment 1
if isfield(tree.subject.frames.frame(1),'position')
    for i=[1:nSamples]
        p_Segment1(i,:)=tree.subject.frames.frame(i).position(1:3);
    end
    %Plot
    figure('name','Position of first segment')
    plot(p_Segment1)
    figure('name','Position of first segment in 3D')
    plot3(p_Segment1(:,1),p_Segment1(:,2),p_Segment1(:,3));
end
