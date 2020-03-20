def generate_inchworm_motor(h_motor, file_name):
    f = open(file_name, "w")
    
    write_prequel()
    
    
    f.close()
    
def write_prequel(f):
    f.write("% angled arm netlist for simulation of statics of angled arm/shuttle interaction")
    f.write("% Developed originally by Daniel S. Contreras in 2018, modified by Ahad Rauf in December 2019")
    
    
    
if name == '__main__':
    generate_inchworm_motor(h_motor, 'subnet_inchworm_motor.m')