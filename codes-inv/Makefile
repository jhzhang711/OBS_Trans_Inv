#-----------------------------------------------------------------------
# Complier  

F77 = mpif90
#-----------------------------------------------------------------------
# Libararies
FFTW   = -I/public1/soft/fftw/3.3.8/include -L/public1/soft/fftw/3.3.8/lib -lfftw3
LAPACK = -llapack -lblas

#-----------------------------------------------------------------------
MOD_OBJS = cls_line_text.o cls_signal_process.o cls_vmodel.o mod_random.o cls_covariance.o cls_recv_func.o cls_disper.o cls_param.o
OBJ_OBJS = dispersion.o lovmrx.o raymrx.o lovdsp.o raydsp.o theo.o qlayer.o four1.o voro2qmodel.o priorvalue.o priorvpvs.o whichcell.o ludcmp.o lubksb.o dsvdcmp.o svbksb.o dpythag.o

# 默认目标
all: driver

# 编译对象文件，确保编译顺序
cls_param.o: cls_param.f90
	$(F77) $(FFLAGS) -c cls_param.f90 -o cls_param.o

cls_line_text.o: cls_line_text.f90
	$(F77) $(FFLAGS) -c cls_line_text.f90 -o cls_line_text.o

cls_signal_process.o: cls_signal_process.f90
	$(F77) $(FFLAGS) ${FFTW} ${LAPACK} -c cls_signal_process.f90 -o cls_signal_process.o


cls_vmodel.o: cls_vmodel.f90 cls_line_text.o
	$(F77) $(FFLAGS) -c cls_vmodel.f90 -o cls_vmodel.o

mod_random.o: mod_random.f90
	$(F77) $(FFLAGS) -c mod_random.f90 -o mod_random.o

cls_covariance.o: cls_covariance.f90
	$(F77) $(FFLAGS) -c cls_covariance.f90 -o cls_covariance.o


###################cls_recv_func.o: cls_vmodel.mod cls_signal_process.mod mod_random.mod cls_covariance.mod
cls_recv_func.o: cls_recv_func.f90 cls_vmodel.o cls_signal_process.o mod_random.o cls_covariance.o
	$(F77) $(FFLAGS) -c cls_recv_func.f90 -o cls_recv_func.o

cls_disper.o: cls_disper.f90 cls_vmodel.o  mod_random.o 
	$(F77) $(FFLAGS) -c cls_disper.f90 -o cls_disper.o

## 其他对象文件的编译规则
%.o: %.f90
	$(F77) $(FFLAGS) $(COV) -c $< $(FFTW) $(LAPACK) -o $*.o 

%.o: %.f
	$(F77) $(FFLAGS) -c $< -o $@


# 生成最终可执行文件
driver: $(MOD_OBJS) $(OBJ_OBJS) Joint.f90
	$(F77) $(FFLAGS) -o JointINV Joint.f90 $(MOD_OBJS) $(OBJ_OBJS) $(FFTW) $(LAPACK) sacio.a -lm

# 清理
clean:
	rm -f *.o *.mod JointINV
