# function zenith(θ, secx, cosx, time)
# end





# subroutine zenith(theta, secx, cosx, ttime)
# double precision:: lat,pi,radian,dec,lha,ttime,theta,cosx,secx
# double precision:: sinld,cosld
# ! solar declination angle from july 1st - harwell traj model
# dec = 23.79
# ! latitude
# lat = 50.0
# pi = 4.0*atan(1.0)
# ! local hour angle - representing time of day
# lha = (1.0+ttime/4.32d+4)*pi
# radian = 180.0/pi
# lat = lat/radian
# dec = dec/radian
# theta = acos(cos(lha)*cos(dec)*cos(lat)+sin(dec)*sin(lat))
# sinld = sin(lat)*sin(dec)
# cosld = cos(lat)*cos(dec)
# cosx = (cos(lha)*cosld)+sinld
# cosx = cos(theta)
# secx = 1.0d+0/(cosx+1.0d-30)
# return

# end

