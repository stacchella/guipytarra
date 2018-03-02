

FC = gfortran
libslatec = /n/sw/fasrcsw/apps/Comp/gcc/4.8.2-fasrc01/slatec/4.1-fasrc01/lib/libslatec.a
libcfitsio = /n/sw/fasrcsw/apps/Comp/gcc/4.8.2-fasrc01/cfitsio/3360-fasrc03/lib/libcfitsio.a
#PGPLOT=
LIBS = $(libcfitsio) $(libslatec)

#ifeq ($(HOSTNAME), orange)
#FC         = gfortran
# PGPLOT =  /usr/lib64/libpgplot.so.5.2.2 /usr/lib64/libpng16.so.16 -lX11 /usr/lib64/libgfortran.so
#libcfitsio = /home/cnaw/lib/libcfitsio.a 
#libslatec  = /home/cnaw/lib/libslatec.a

#LIBS       =  $(libcfitsio) $(libslatec) -L/usr/local/lib64/libquadmath.so.0
#endif

FFLAGS = -O3 -g -pedantic -C -mcmodel=medium 
#FFLAGS = -O3 -g  -L/usr/local/lib/libmpc.so.3

AR  = ar
RM  = rm -f
SOURCES = ab_mag_to_photon_flux.f \
	acs_cosmic_rays.f \
	add_baseline.f \
	add_dark.f \
	add_galaxy_component.f \
	add_ipc.f \
	add_latents.f \
	add_modelled_cosmic_rays.f \
	add_modelled_galaxy.f \
	add_one_over_f_noise.f \
	add_read_noise.f \
	add_reference_pixels.f \
	add_sky_background.f \
	add_stars.f \
	add_up_the_ramp.f \
	alngam.f \
	clear_accum.f \
	closefits.f \
	coadd.f \
	cr_distribution.f \
	data_model_fits.f \
	divide_frame.f \
	find_sca_coords.f \
	gain_per_amp.f \
	gammain_bhattacharjee.f \
	get_date.f \
	getkeywords.f \
	getfitsdata2d.f \
	get_sca_id.f \
	integration.f \
	inverse_correction.f \
	jl_latent.f \
	julian_day.f \
	jwst_keywords.f \
	linearity_incorrect.f \
	linear_interpolation.f \
	linearity_table.f \
	load_osim_transforms.f \
	make_fake_cat.f \
	median.f \
	open_big_fits_cube.f \
	openfits.f \
	osim_coords_from_sca.f \
	plane_to_spherical.f \
	printerror.f \
	psf_convolve.f \
	randgen.f \
	ra_dec_to_sca.f \
	ra_dec_xy_osim.f \
	read_cr_matrix.f \
	read_fake_mag_cat.f \
	read_filter_parameters.f \
	read_fits.f \
	read_funky_fits.f \
	read_jwst_background.f \
	read_psf.f \
	read_linearity.f \
	read_psf_list.f \
	read_psf_fits.f \
	read_sca_calib.f \
	read_sca_wcs.f \
	read_sim_cat.f \
	read_star_catalogue.f \
	resample.f \
	rotate_coords.f \
	sca_boundaries.f \
	sca_coords_from_osim.f \
	sca_footprint.f \
	sca_to_ra_dec.f \
	set_params.f \
	smooth.f \
	spherical_to_plane.f \
	spline_interp.f \
	structural.f \
	total_time.f \
	ut_to_date_time.f \
	wcs_keywords.f \
	wcs_keywords_new.f \
	write_2d.f \
	write_float_2d_image.f \
	write_float_3d_image.f \
	write_float_image_section.f \
	write_int_2d_image.f \
	write_frame.f \
	write_int_3d_image.f \
	write_int_image_section.f \
	write_nircam_keywords.f \
	write_uint_image_section.f \
	xbiwt.f \
	xy_osim_ra_dec.f \
	zeroin.f \
	zodi.f \
	zodi_background.f 
#
OBJECTS := $(SOURCES:.f=.o)

capenga: capenga.o $(OBJECTS)
	$(FC) $(FFLAGS) -o capenga \
	capenga.o $(OBJECTS) $(LIBS)
#
create_star_catalogue: create_star_catalogue.o $(OBJECTS)
	$(FC) $(FFLAGS) -o create_star_catalogue \
	create_star_catalogue.o $(OBJECTS) \
	$(LIBS)
#
debris_disk: debris_disk.o frebin.o read_fits.o openfits.o closefits.o \
printerror.o write_fits_image.o getfitsdata2d.o integration.o \
linear_interpolation.o randgen.o read_psf.o read_psf_fits.o psf_convolve.o \
resample.o add_ipc.o
	$(FC)  $(FFLAGS) -o debris_disk debris_disk.o \
	frebin.o read_fits.o openfits.o closefits.o \
	printerror.o write_fits_image.o getfitsdata2d.o integration.o \
	linear_interpolation.o randgen.o \
	read_psf.o read_psf_fits.o psf_convolve.o \
	resample.o add_ipc.o \
	$(LIBS)

find_zp: find_zp.o read_filter_parameters.o ab_mag_to_photon_flux.o \
zeroin.o
	$(FC) $(FFLAGS) -o find_zp find_zp.o \
	read_filter_parameters.o ab_mag_to_photon_flux.o \
	zeroin.o
#
guitarra: guitarra.o $(OBJECTS) ; echo $(FC)
	$(FC) $(FFLAGS) -o guitarra \
	guitarra.o $(OBJECTS) $(LIBS)
#
#match: match.o catmatch.o read_cat.o flgrid.o vector_plot.o mystat.o
#	$(FC) $(FLAGS) -o match match.o \
#	catmatch.o read_cat.o flgrid.o vector_plot.o mystat.o \
#	$(PGPLOT) -lX11 $(LIBPNG)
#
proselytism: proselytism.o ra_dec_to_sca.o load_osim_transforms.o \
ra_dec_xy_osim.o sca_coords_from_osim.o rotate_coords.o spherical_to_plane.o \
read_nircam_outline.o rotate_coords2.o
	$(FC) $(FFLAGS) -o proselytism proselytism.o \
	ra_dec_to_sca.o load_osim_transforms.o \
	ra_dec_xy_osim.o sca_coords_from_osim.o rotate_coords.o spherical_to_plane.o \
	read_nircam_outline.o rotate_coords2.o
#
psf_loss: psf_loss.o write_test_image.o $(OBJECTS)
	$(FC) $(FFLAGS) -o psf_loss \
	psf_loss.o write_test_image.o $(OBJECTS) $(LIBS)
#
simulator : simulator.o $(OBJECTS) ; echo $(FC)
	$(FC) $(FFLAGS) -o simulator \
	simulator.o $(OBJECTS) $(LIBS)
#
#simulate_phot_cat : simulate_phot_cat.o $(OBJECTS) 
#	$(FC) $(FFLAGS) -o simulate_phot_cat \
#	plot_sed.o plot_set_colors.o \
#	simulate_phot_cat.o $(OBJECTS) $(LIBS) $(PGPLOT)
#
test_data_model : test_data_model.o \
data_model_fits.o jwst_keywords.o $(OBJECTS)
	$(FC) $(FFLAGS) -o test_data_model \
	test_data_model.o data_model_fits.o jwst_keywords.o \
	$(OBJECTS) $(LIBS)
# 
test_make_fake_cat : test_make_fake_cat.o \
make_fake_cat.o \
read_fake_mag_cat.o \
read_candels_goods_n.o \
read_candels_goods_s.o
#
	$(FC) $(FFLAGS) -o test_make_fake_cat \
	test_make_fake_cat.o \
	make_fake_cat.o \
	read_fake_mag_cat.o \
	read_candels_goods_n.o \
	read_candels_goods_s.o \
	$(LIBS)
#
test_psf : test_psf.o write_test_image.o $(OBJECTS)
	$(FC) $(FFLAGS) -o test_psf test_psf.o \
	write_test_image.o \
	$(OBJECTS) $(LIBS)
#
test_randgen : test_randgen.o randgen.o median.o xbiwt.o
	$(FC) $(FFLAGS) -o test_randgen test_randgen.o \
	randgen.o median.o xbiwt.o $(libslatec)

test_total_time: test_total_time.o total_time.o set_params.o
	$(FC) $(FFLAGS) -o test_total_time test_total_time.o \
	total_time.o set_params.o
#
test_wcs_keywords: test_wcs_keywords.o $(OBJECTS)
	$(FC) $(FFLAGS) -o test_wcs_keywords test_wcs_keywords.o \
	$(OBJECTS) $(LIBS)
#
.PHONY : clean
clean : 
	$(RM) $(OBJECTS)
