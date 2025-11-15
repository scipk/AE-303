%% Description -- Please Read!
%
% This script imports the data from the .csv files for the PIV lab of
% AE303. A status is printed to the command window as each file is imported
% and then once again as each file's data is formatted. Each .csv file
% represents one instant in time for the flow. The data was recorded at 150
% Hz, so 300 files corresponds to 2 seconds of data. 
%
% This script checks if the data has been read before so that one does not
% spend extraneous time re-importing data. However, if the import portion
% of the script is terminated early, the parsing portion of the script will
% not function properly and the workspace should be cleared before running
% this script again. 
%
% To further avoid repeating portions of the script, divide the script
% into sections with the '%% Title' as is done already and take advantage
% of the 'Run Section' option to the right of 'Run' in Matlab.
%
% The purpose of this script is to allow the student to work with the data
% presented rather than put unnecessary time into writing their own import
% and formatting script since there is such a large amount of data.
%
%% Notes on the script output
%
% Upon this script's completion, there will exist 10 new variables in the
% workspace. Each is important for processing the PIV data. They are as
% follows:
%
%   data        This is the raw 3D array containing all of the data from
%               all 300 of the .csv files; dimensions 8806x11x300
%
%   N           The number of time steps in the dataset; 300
%
%   winsize     The 1x2 array with the number of rows and columns,
%               respectively, in the vector grid; should be [119,74]
%
%   X           The 2D matrix of the x-coordinates of the vector grid
%
%   Y           The 2D matris of the y-coordinates of the vector grid
%
%   u           The 2D matrix of the u-component of the velocity on the
%               vector grid with likely erroneous vectors replaced by NaN
%               with units of m/s
%
%   v           The 2D matrix of the v-component of the velocity on the
%               vector grid with likely erroneous vectors replaced by NaN
%               with units of m/s
%
%   u_raw       The 2D matrix of the u-component of the velocity on the
%               vector grid with no replaced vectors and units of m/s
%
%   v_raw       The 2D matrix of the v-component of the velocity on the
%               vector grid with no replaced vectors and units of m/s
%
%   u_pix       The 2D matrix of the u-component of the velocity on the
%               vector grid with likely erroneous vectors replaced by NaN
%               with units of pixel displacement
%
%   v_pix       The 2D matrix of the v-component of the velocity on the
%               vector grid with likely erroneous vectors replaced by NaN
%               with units of pixel displacement
%
% The velocity output arrays are all 3D arrays, where the first two
% dimensions represent the vector grid and the third dimension is how each
% time step is stored. Moving, for example, from u(:,:,1) to u(:,:,2) is 
% moving one time step forward in time, while referencing the full vector
% grid of u velocity components.
%
% The bottom of the parse portion of the script includes a commented out
% code which will plot and animation of the flow as the data is parsed.
% Running this portion can cause the code to take longer than desired to
% run, which is why it is commented out.
%
%% Read Data
winsize = [119 74]; % Size of vector grid
L = winsize(1)*winsize(2); % Get size of vector field

if ~exist('data','var') % Check if data read already

    N = 300; % Number of files (time steps)
    
    data = zeros(L,11,N); % Pre-allocate

    for i = 1:N
        fprintf('Now Retrieving File %d of 300\n',i);
        file = sprintf('EduPIV_lab.62tbxosb.%06d.csv',i-1);
        data(:,:,i) = readmatrix(file); % Save each file to data
    end
end

%% Parse Data
s = 1;
t = 1;
u = zeros(119,74,300); v = u;
u_raw = u; u_pix = u;
v_raw = v; v_pix = v; % Pre-Allocate
x = u; y = v;
for i = 1:N
    fprintf('Now Parsing Set %d of 300\n',i); % Status display
    for q = 1:8806
        if data(q,11,i) ~= 0
            u(s,t,i) = nan; % Replace flagged vectors with nan
            v(s,t,i) = nan;
            u_pix(s,t,i) = nan;
            v_pix(s,t,i) = nan;
        else
            u(s,t,i) = data(q,9,i); % Parse valid vectors
            v(s,t,i) = data(q,10,i);
            u_pix(s,t,i) = data(q,7,i);
            v_pix(s,t,i) = data(q,8,i);
        end
        
        u_raw(s,t,i) = data(q,9,i);
        v_raw(s,t,i) = data(q,10,i);
        
        x(s,t,i) = data(q,5,i); % Record grid
        y(s,t,i) = data(q,6,i);
        
        s = s + 1; % Move to next row
        
        if mod(q,winsize(1)) == 0 % Detect is done with column
            t = t + 1; % Move to next vector column
            s = 1; % Reset to first row
        end
        if q == L
            t = 1; % Finished with frame, reset t
            % figure(1) % Update only figure 1
            % contourf(x(:,:,i),y(:,:,i),sqrt(u(:,:,i).^2 + v(:,:,i).^2),20,'LineStyle','none')
            % axis equal % Update cool animation of the flow
            % xlabel('x [cm]')
            % ylabel('y [cm]')
            % drawnow
        end
    end
end

X = x(:,:,1); % Grid doesn't change, so we just need a 2D matrix
Y = y(:,:,1);

clear file i q s t x y L

%% Data Processing





