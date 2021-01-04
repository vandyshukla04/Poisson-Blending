%Vandita Shukla 17141590
clear all;

blend_this_image = double(rgb2gray(imread('images/womanyellingcat.jpg')))/255 ;

%we find the size of the image
[R,C] = size(blend_this_image);

region_to_blend = roipoly(blend_this_image);
region_to_blend = double(region_to_blend);

image1 = blend_this_image.*(double(ones(R,C)-region_to_blend));

imshow(image1);

%we will first count the pixels in mask
mask_pixel_count = 0;
for i=1:R
    for j=1:C
        if(region_to_blend(i,j)==1)
            region_to_blend(i,j)=region_to_blend(i,j)+mask_pixel_count;
            mask_pixel_count = mask_pixel_count+1;
        end
    end
end

omega_number = mask_pixel_count;

%since we are using a 4-neighbour model and we also need to test the
%boundaries of the image, we will create the boundary of Omega i.e. the cropped out
%region after padding the region_to_blend

region_to_blend_2 = padarray(region_to_blend,[1,1]);
[M,N] = size(region_to_blend_2);
boundary = zeros(R,C);

isboundary=0;

for i=2:M-1
    for j=2:N-1
        %3x3 patch for testing 4 neighbours
        patch = region_to_blend_2(i-1:i+1,j-1:j+1);
        %test whether center of patch is a boundary pixel
        if(patch(2,2)==0)
            if(patch(2,1)>0)
                isboundary=1;
            elseif(patch(2,3)>0)
                isboundary=1;
            elseif(patch(1,2)>0)
                isboundary=1;
            elseif(patch(3,2)>0)
                isboundary=1;
            end
        end
        %if it is a boundary patch then fill in the boundary array
        if(isboundary==1)
            boundary(i-1,j-1)=1;
        end
        %reset the flag
        isboundary=0;      
    end
end

%now we get to the task of filling the b array
%we pad boundary array and the region_to_blend with 1 extra row and column on every border since we
%will use 4 neighbour model
%b array contains in this task no v_pq because guidance of g is 0 in this
%particular case. So b only is the sum of all neighbouring pixels that are
%a part of the boundary.

boundary = padarray(boundary,[1,1]);
region_to_blend = padarray(region_to_blend,[1,1]);

b = zeros(omega_number,1);

for i = 2:M-1
    for j = 2:N-1
        
        if(region_to_blend(i,j)>0)
            %let the numbered pixels be the counter and sum be for
            %intensities
            counter=region_to_blend(i,j);
            sum=0;
            
            if(boundary(i-1,j))
                sum = sum + blend_this_image(i-2,j-1);
            end
            if(boundary(i+1,j))
                sum = sum + blend_this_image(i,j-1);
            end
            if(boundary(i,j-1))
                sum = sum + blend_this_image(i-1,j-2);
            end
            if(boundary(i,j+1))
                sum = sum + blend_this_image(i-1,j-2);
            end
            
            b(counter) = sum;
        end
    end
end

%Now we create the A matrix and calculate the diagonal
%For calculating A matrix we will keep a counter for a pixel from mask and
%all it's four neighbours. A will also be a sparse matrix. This will require us to create arrays
%5 times the size of omega_number i.e. the number of pixels in the region
%omega. 

s_i = zeros(omega_number*5,1);%storing the value of pixel whose neighbours we are calculating
s_j = zeros(omega_number*5,1);%storing which are the neighbours of the current pixel in the region
%excluding boundary
v = zeros(omega_number*5,1);%will store -1 for existent neighbour pixels within omega to calculate laplacian
position = 1;

for i = 2:M-1
    for j = 2:N-1
        
        if(region_to_blend(i,j)>0)
            counter=region_to_blend(i,j);
        
            s_i(position,1)=counter;
            s_j(position,1)=counter;
            v(position,1)=4;%Np diagonal element
            
            %now we fill the rest array
            position = position+1;
            
            if(region_to_blend(i-1,j)>0)
                s_i(position,1)=counter;
                s_j(position,1)=region_to_blend(i-1,j);
                v(position,1)=-1;
                position=position+1;
            end
            
            if(region_to_blend(i+1,j)>0)
                s_i(position,1)=counter;
                s_j(position,1)=region_to_blend(i+1,j);
                v(position,1)=-1;
                position=position+1;
            end
            
            if(region_to_blend(i,j-1)>0)
                s_i(position,1)=counter;
                s_j(position,1)=region_to_blend(i,j-1);
                v(position,1)=-1;
                position=position+1;
            end
            
            if(region_to_blend(i,j+1)>0)
                s_i(position,1)=counter;
                s_j(position,1)=region_to_blend(i,j+1);
                v(position,1)=-1;
                position=position+1;
            end
            
        end
    end
end

%we now reduce array till last filled position and create sparse matrix A

s_i=s_i(1:position-1,:);
s_j=s_j(1:position-1,:);
v=v(1:position-1,:);

A = sparse(s_i,s_j,v,omega_number,omega_number);

x=A\b;

%now filling the empty region in image1

for i=2:M-2
    for j=2:N-2
        if(region_to_blend(i,j)>0)
            image1(i-1,j-1)=x(region_to_blend(i,j),1);
        end
    end
end

imshow(image1);
saveas(gcf, 'task_1.jpg') ;
            