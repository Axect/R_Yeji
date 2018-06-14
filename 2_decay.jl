# The rate of the two-body decay channels

using Plots

# mass of particle

m_H = 125.09
m_Z = 91.1876
m_W = 80.385
m_fu = 2.2 * 1e-3
m_fd = 4.7 * 1e-3
m_fs = 96 * 1e-3
m_fc = 1.28
m_ft = 173.1
m_fb = 4.18
m_e = 0.5109989461 * 1e-3
m_mu = 105.6583745 * 1e-3
m_tau = 1776.86 * 1e-3
m_ne = 1.85 * 1e-9
m_nmu = 0.19 * 1e-3
m_ntau = 18.2 * 1e-3

# the number of color charge

N_c = 3

function genYH(m_DM :: Float64)
	if m_DM <= 2*m_H
		return 0.
	else
  	x_H = (m_H/ m_DM)^2
  	return m_DM^(3)*(1 + 2*x_H)^(2)*(1 - 4*x_H)^(1/2)
  end
end

function genYZ(m_DM :: Float64)
	if m_DM <= 2*m_Z
		return 0.
	else
	  x_Z = (m_Z/ m_DM)^2
  	return m_DM^(3)*(1 - 4x_Z + 12*x_Z^(2))*(1 - 4*x_Z)^(1/2)
	end
end

function genYW(m_DM :: Float64)
	if m_DM <= 2*m_W
		return 0.
	else
  	x_W = (m_W/ m_DM)^2
  	return 2* m_DM^(3)*(1 - 4*x_W + 12*x_W^(2))*(1 - 4*x_W)^(1/2)
	end
end

function genYf(m_f, m_DM :: Float64)
  x_f = (m_f/ m_DM).^2
  return N_c * 4* m_DM * 4* m_f^2 * (1-4*x_f)^(3/2)
end


yH = map(x -> genYH(x), 1:0.25:1e6)
yZ = map(x -> genYZ(x), 1:0.25:1e6)
yW = map(x -> genYW(x), 1:0.25:1e6)
# yfu = map(x -> genYf(m_fu, x), filter(y -> y > m_fu*2, 1:0.25:1e6))
# yfd = map(x -> genYf(m_fd, x), filter(y -> y > m_fd*2, 1:0.25:1e6))
# yfs = map(x -> genYf(m_fs, x), filter(y -> y > m_fs*2, 1:0.25:1e6))
# yfc = map(x -> genYf(m_fc, x), filter(y -> y > m_fc*2, 1:0.25:1e6))
# yft = map(x -> genYf(m_ft, x), filter(y -> y > m_ft*2, 1:0.25:1e6))
# yfb = map(x -> genYf(m_fb, x), filter(y -> y > m_fb*2, 1:0.25:1e6))
# ye = map(x -> genYf(m_e, x), filter(y -> y > m_e*2, 1:0.25:1e6))
# ymu = map(x -> genYf(m_mu, x), filter(y -> y > m_mu*2, 1:0.25:1e6))
# ytau = map(x -> genYf(m_tau, x), filter(y -> y > m_tau*2, 1:0.25:1e6))
# yne = map(x -> genYf(m_ne, x), filter(y -> y > m_ne*2, 1:0.25:1e6))
# ynmu = map(x -> genYf(m_nmu, x), filter(y -> y > m_nmu*2, 1:0.25:1e6))
# yntau = map(x -> genYf(m_ntau, x), filter(y -> y > m_ntau*2, 1:0.25:1e6))


#function Yf()
#	yy = collect(1:0.25:1e6)
#	for (i, y) in enumerate(1:0.25:1e6)
#		if y < 2*m_ne
#			yy[i] = 0
#		elseif y < 2*m_nmu
#			yy[i] = genYf(m_ne, y)
#		elseif y < 2*m_e
#			yy[i] = genYf(m_nmu, y) + genYf(m_ne, y)
#		elseif y < 2*m_fu
#			yy[i] = genYf(m_e, y) + genYf(m_nmu, y) + genYf(m_ne, y)
#		elseif y < 2*m_fd	
#			yy[i] = genYf(m_fu, y) + genYf(m_e, y) + genYf(m_nmu, y) + genYf(m_ne, y)
#		elseif y < 2*m_ntau
#			yy[i] = genYf(m_fd, y) + genYf(m_fu, y) + genYf(m_e, y) + genYf(m_nmu, y) + genYf(m_ne, y)
#		elseif y < 2*m_fs
#			yy[i] = genYf(m_ntau, y) + genYf(m_fd, y) + genYf(m_fu, y) + genYf(m_e, y) + genYf(m_nmu, y) + genYf(m_ne, y)
#		elseif y < 2*m_mu
#			yy[i] = genYf(m_fs, y) + genYf(m_ntau, y) + genYf(m_fd, y) + genYf(m_fu, y) + genYf(m_e, y) + genYf(m_nmu, y) + genYf(m_ne, y)
#		elseif y < 2*m_fc
#			yy[i] = genYf(m_mu, y) + genYf(m_fs, y) + genYf(m_ntau, y) + genYf(m_fd, y) + genYf(m_fu, y) + genYf(m_e, y) + genYf(m_nmu, y) + genYf(m_ne, y)
#		elseif y < 2*m_tau	
#			yy[i] = genYf(m_fc, y) + genYf(m_mu, y) + genYf(m_fs, y) + genYf(m_ntau, y) + genYf(m_fd, y) + genYf(m_fu, y) + genYf(m_e, y) + genYf(m_nmu, y) + genYf(m_ne, y)
#		elseif y < 2*m_fb
#			yy[i] = genYf(m_tau, y)+ genYf(m_fc, y) + genYf(m_mu, y) + genYf(m_fs, y) + genYf(m_ntau, y) + genYf(m_fd, y) + genYf(m_fu, y) + genYf(m_e, y) + genYf(m_nmu, y) + genYf(m_ne, y)
#		elseif y < 2*m_ft
#			yy[i] = genYf(m_fb, y) + genYf(m_tau, y) + genYf(m_fc, y) + genYf(m_mu, y) + genYf(m_fs, y) + genYf(m_ntau, y) + genYf(m_fd, y) + genYf(m_fu, y) + genYf(m_e, y) + genYf(m_nmu, y) + genYf(m_ne, y)
#		else
#			yy[i] = genYf(m_ft, y) + genYf(m_fb, y) + genYf(m_tau, y) + genYf(m_fc, y) + genYf(m_mu, y) + genYf(m_fs, y) + genYf(m_ntau, y) + genYf(m_fd, y) + genYf(m_fu, y) + genYf(m_e, y) + genYf(m_nmu, y) + genYf(m_ne, y)
#		end
#	end
#	return yy
#end

#yff = Yf()



#body2= yH + yZ + yW + yff 

plotly();

plot(1:0.25:1e06, yH, xaxis=:log10, yaxis=:log10)
#plot(1:0.25:1e06, yZ+yW, xaxis=:log10, yaxis=:log10)
#plot(1:0.25:1e06, yW, xaxis=:log10, yaxis=:log10)
#plot(1:0.25:1e06, yff, xaxis=:log10, yaxis=:log10)

