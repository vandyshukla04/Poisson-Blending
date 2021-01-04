%Vandita Shukla 17141590
clear all;

%blend_this_image = double(rgb2gray(imread('images/dsc_3656.jpg')))/255 ;
%blend_this_image = imresize(blend_this_image,0.4);
%get_region_from = double(rgb2gray(imread('images/animals1.jpg')))/255 ;

blend_this_image = double(rgb2gray(imread('images/mars.jpg')))/255 ;
get_region_from = double(rgb2gray(imread('images/lion.jpg')))/255 ;
get_region_from = imresize(get_region_from,0.9);

[RB, CB] = size(blend_this_image);
[RG, CG] = size(get_region_from);
%we find the size of the image

[mask1, x, y] = roipoly(get_region_from);
mask1 = double(mask1);
%image1 = get_from.*(double(ones(R,C)-mask1));

%now we select the origin pixel of the region map we obtained from blend_this_image 
[col,row,~]= impixel(blend_this_image);
x = x - min(x);
y = y - min(y);
x = x + col;
y = y + row;

mask2 = roipoly(blend_this_image,x,y);
imshow(mask2);

mask2 = double(mask2);
image1 = blend_this_image.*(double(ones(RB,CB)-mask2));
imshow(image1);

%use this mask
%we will first count the pixels in mask
mask_pixel_count = 0;
for i=1:RB
    for j=1:CB
        if(mask2(i,j)==1)
            mask2(i,j)=mask2(i,j)+mask_pixel_count;
            mask_pixel_count = mask_pixel_count+1;
        end
    end
end

omega_number = mask_pixel_count;

%since we are using a 4-neighbour model and we also need to test the
%boundaries of the image, we will create the boundary of Omega i.e. the cropped out
%region after padding the region_to_blend

region_to_blend_2 = padarray(mask2,[1,1]);
[M,N] = size(region_to_blend_2);
boundary = zeros(RG,CG);

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

%now we get to the task of filling the b array and this time b array will
%also have vpq

laplacian = [0 -1 0;-1 4 -1;0 -1 0];
convoluted = conv2(get_region_from,laplacian,'same');

vpq = zeros(1,RG*CG);

vpq_counter = 1;

for i = 1:RG
    for j = 1:CG
        if(mask1(i,j)>0)
            vpq(1,vpq_counter)= convoluted(i,j);
            vpq_counter = vpq_counter+1;
        end
    end
end

%we pad boundary array and the region_to_blend with 1 extra row and column on every border since we
%will use 4 neighbour model

boundary = padarray(boundary,[1,1]);
mask2 = padarray(mask2,[1,1]);

b = zeros(omega_number,1);

for i = 2:M-1
    for j = 2:N-1
        
        if(mask2(i,j)>0)
            %let the numbered pixels be the counter and sum be for
            %intensities
            counter=mask2(i,j);
            sum=vpq(1,counter);
            
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
        
        if(mask2(i,j)>0)
            counter=mask2(i,j);
        
            s_i(position,1)=counter;
            s_j(position,1)=counter;
            v(position,1)=4;%Np diagonal element
            
            %now we fill the rest array
            position = position+1;
            
            if(mask2(i-1,j)>0)
                s_i(position,1)=counter;
                s_j(position,1)=mask2(i-1,j);
                v(position,1)=-1;
                position=position+1;
            end
            
            if(mask2(i+1,j)>0)
                s_i(position,1)=counter;
                s_j(position,1)=mask2(i+1,j);
                v(position,1)=-1;
                position=position+1;
            end
            
            if(mask2(i,j-1)>0)
                s_i(position,1)=counter;
                s_j(position,1)=mask2(i,j-1);
                v(position,1)=-1;
                position=position+1;
            end
            
            if(mask2(i,j+1)>0)
                s_i(position,1)=counter;
                s_j(position,1)=mask2(i,j+1);
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
        if(mask2(i,j)>0)
            image1(i-1,j-1)=x(mask2(i,j),1);
        end
    end
end

imshow(image1);
saveas(gcf, 'task_2a.jpg') ;
            