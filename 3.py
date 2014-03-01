import numpy as np
import pyfits as pf
import matplotlib.pylab as plt

file_name='/Users/anita/Documents/University_Third_Year/AST326/Lab5/GJ1214/star'
centroid_name='/Users/anita/Documents/University_Third_Year/AST326/Lab5/all_centroid'
star_list=np.arange(0,364,1)
#index=[9]#,20,38,82,84,88,103,104,133,145,149,161,171,250,276,282,359]
#star_list=np.delete(star_list_temp,index)


bias1='/Users/anita/Documents/University_Third_Year/AST326/Lab5/bias.fit'
bias_array_final=pf.getdata(bias1)
flat1='/Users/anita/Documents/University_Third_Year/AST326/Lab5/flat.fit'
flat_array=pf.getdata(flat1)
S_N=[]
Noise=[]
fflux=[]
bbackground=[]
averaged_list=[]
normalized_list=[]

for j in star_list:
    stars=file_name+str(j)+'.fit'
    star=pf.getdata(stars)
    corrected=(star-bias_array_final)/((flat_array))
    temp_flux=np.zeros(35)
    temp_background=[]
    temp_flux=[]
    total_noise_list=[]
    signal_to_noise_list=[]
    centroid_x_list=np.loadtxt('/Users/anita/Documents/University_Third_Year/AST326/Lab5/all_centroid/'+str(j)+'centroids.txt',usecols=(0,))
    centroid_y_list=np.loadtxt('/Users/anita/Documents/University_Third_Year/AST326/Lab5/all_centroid/'+str(j)+'centroids.txt',usecols=(1,))
    #print "Shape for star number:",j,"is:",(np.shape(centroid_x_list))
    #output=np.zeros((1,80))
    print "star number:", j
    i =0
    while i < 20:
        r =centroid_x_list[i]
        u =centroid_y_list[i]
    	box_x=np.arange(-26,26,1)  #box in x direction
    	box_y=np.arange(-26,26,1)  #box in y direction
    	total_flux=[]
    	flux_list=[]
    	background_list=[]
    	for l in (box_x+r):
            for p in (box_y+u):
                distance=np.sqrt(((l-r)**2)+((p-u)**2)) #distance between centroid and the points in the box
                print "distance is:", distance
                if distance <15:
                    flux_list.append(corrected[l][p]) #appends the flux of the point if inside radius 15
            	elif distance<25 and distance>22:
                    background_list.append(corrected[l][p])  #appends the background if inside annulus(22<r<25)
    	flux_sum=np.sum(flux_list)
    	temp_flux.append(flux_sum)
    	background_avg=np.mean(background_list)
    	temp_background.append(background_avg)
    	i+=1
    # noise
    output=np.column_stack((temp_flux,temp_background))
    np.savetxt('fb{0}.txt'.format(j),output,fmt='%.1i')
    fflux.append(temp_flux)
    bbackground.append(temp_background)
    background_array=np.array(bbackground)
    r1=15
    r2=22
    r3=25
    readout_noise=6.0/1.9
    N1=(np.pi)*(r1**2)
    N23=(np.pi)*((r3**2)-(r2**2))
    B=background_array/(np.size(bbackground))
    poisson_noise=np.std(fflux)
    total_noise=(poisson_noise)+(N1*(B+(readout_noise**2)))+((N1*(B+(readout_noise**2)))/N23)
    total_noise_list.append(total_noise)
    signal_to_noise=poisson_noise/(np.sqrt(total_noise))
    signal_to_noise_list.append(signal_to_noise)
    w=1./total_noise[j]
    numerator=np.sum(w*fflux[j])
    denominator=np.sum(w)
    averaged=numerator/denominator
    averaged_list.append(averaged)
    normalized=fflux[j]/averaged_list[j]
    normalized_list.append(normalized)
    #output1=np.column_stack((total_noise,signal_to_noise))
    #np.savetxt('properties{0}.txt'.format(j),output,fmt='%.1i')

single_flux=[]
double_flux=[]
triple_flux=[]
quadra_flux=[]
qwe_list=[]

for j in star_list:
    stars=file_name+str(j)+'.fit'
    star=pf.getdata(stars)
    corrected=(star-bias_array_final)/((flat_array))
    ffname='/Users/anita/Documents/University_Third_Year/AST326/Lab5/fb'+str(j)+'.txt'
    flux_file=np.loadtxt(ffname,usecols=(0,))
    flux_file=flux_file[::-1]
    #single_flux.append(flux_file[0])
    #double_flux.append(flux_file[1])
    #triple_flux.append(flux_file[2])
    #quadra_flux.append(flux_file[3])
    qwe=normalized_list[j][0]
    qwe_list.append(qwe)
                       
time_array=np.linspace(2455704.4993055,2455704.6843287,364)
#plt.subplot(4,5,m)
plt.plot(time_array,qwe_list,'b.')
#plt.plot(time_array,double_flux,'m.')
#plt.plot(time_array,triple_flux,'r.')
#plt.plot(time_array,quadra_flux,'g.')
plt.show()




    