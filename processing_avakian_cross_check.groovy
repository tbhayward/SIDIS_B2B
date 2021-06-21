/*
 * author Timothy B. Hayward
 * 
 * SIDIS dihadron 
 */

import java.io.File;

import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;
import org.jlab.clas.physics.*;
import org.jlab.clas12.physics.*;

// import from hayward_coatjava_extensions
import extended_kinematic_fitters.*; 
import analyzers.*;

// dilks CLAS QA analysis
import clasqa.QADB

// filetype for gathering files in directory
import groovy.io.FileType;


public class processing_avakian_cross_check {

	public static double phi_calculation (double x, double y) {
		// tracks are given with Cartesian values and so must be converted to cylindrical
		double phi = Math.toDegrees(Math.atan2(x,y));
		phi = phi - 90;
		if (phi < 0) {
			phi = 360 + phi;
		}
		phi = 360 - phi;
		return phi;	
	}

	public static void main(String[] args) {
		File[] hipo_list;
		if (args.length == 0) {
			// exits program if input directory not specified 
    	   	println("ERROR: Please enter a hipo file directory as the first argument");
    	  	System.exit(0);
    	} else {
    		File directory = new File(args[0]);
    		hipo_list = directory.listFiles();
    	}

    	def list = []

		def dir = new File(args[0])
		dir.eachFileRecurse (FileType.FILES) { file ->
		  list << file
		  // println(file.toString()); println(); 
		}

		println(); println(); println();

		String p1_Str;
		if (args.length < 2) {
			// assigns pi+ to p1
			println("WARNING: Specify a PDG PID for p1! Set to pi+. \n");
			p1_Str = "211";
		} else {
			p1_Str = args[1];
			println("Set p1 PID = "+p1_Str+"\n");
		}

		String p2_Str;
		if (args.length < 3) {
			// assigns pi- to p2
			println("WARNING: Specify a PDG PID for p2! Set to pi-. \n");
			p2_Str = "-211";
		} else {
			p2_Str = args[2];
			println("Set p2 PID = "+p2_Str+"\n");
		}

		String output_file;
		if (args.length < 4) {
			// uses dummy name for output file if not specified
			println("WARNING: Specify an output file name. Set to \"dihadron_dummy_out.txt\".\n");
			output_file = "dihadron_dummy_out.txt"
		} else {
			output_file = args[3];
		}

		int n_files;
		if ((args.length < 5)||(Integer.parseInt(args[4])>hipo_list.size())) {
			// if number of files not specified or too large, set to number of files in directory
			println("WARNING: Number of files not specified or number too large."); 
			println("Setting # of files to be equal to number of files in the directory.");
			// n_files = hipo_list.size();
			n_files = list.size();
			println("There are "+hipo_list.size()+" or maybe "+list.size()+" number of files.")
		} else{
			// if specified, convert to int
			n_files = Integer.parseInt(args[4]);
		}

		File file = new File(output_file);
		file.bytes = new byte[0]

		int hadron_pair_counts = 0;
		GenericKinematicFitter research_fitter = new analysis_fitter(10.6041); 
		// load my kinematic fitter/PID
		EventFilter filter = new EventFilter("11:"+p1_Str+":"+p2_Str+":X+:X-:Xn"); 
		// set filter for final states
		// setup QA database
		QADB qa = new QADB();

		int num_events = 0;
		int current_file = 0;
		// for (int current_file; current_file<n_files; current_file++) {
		while (current_file < n_files) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files); println(); println();
			// limit to a certain number of files defined by n_files

			HipoDataSource reader = new HipoDataSource();

			reader.open(list[current_file]); // open next hipo file
			current_file++;
			HipoDataEvent event = reader.getNextEvent(); 

			while(reader.hasEvent()==true){
			// for (int i=0; i<100000; i++) {
				num_events++; 
				if (num_events%50000 == 0) { // not necessary
					print("processed: "+num_events+" events. ");
				}

				// get run and event numbers
				event = reader.getNextEvent();
			    int runnum = event.getBank("RUN::config").getInt('run',0); // collect info for QA
			    int evnum = event.getBank("RUN::config").getInt('event',0);

			    PhysicsEvent research_Event = research_fitter.getPhysicsEvent(event);

			    boolean process_event = false;
			    // if (runnum == 11) { // if run number = 11 then it is MC and we don't use QA
			    // 	process_event = filter.isValid(research_Event);
			    // } else {
			    // 	process_event = (filter.isValid(research_Event) && 
			    // 		qa.OkForAsymmetry(runnum,evnum));
			    // }
			    // going to do QA as a status variable
			    process_event = filter.isValid(research_Event);
				if (process_event) {
					int num_p1 = research_Event.countByPid(p1_Str.toInteger());  
					// get # of particles w/ pid1
					int num_p2 = research_Event.countByPid(p2_Str.toInteger()); 
					// get # of particles w/ pid2

					// num_p1 = 1; num_p2 = 1;
					for (int current_p1 = 0; current_p1 < num_p1; current_p1++) { 
					// cycle over all combinations
						for (int current_p2 = 0; current_p2 < num_p2; current_p2++) {
					// for (int current_p1 = 0; current_p1 < 1; current_p1++) { 
					// // cycle over all combinations
					// 	for (int current_p2 = 0; current_p2 < 1; current_p2++) {

							Dihadrons variables = new Dihadrons(event, research_Event, 
								p1_Str.toInteger(), current_p1, p2_Str.toInteger(), current_p2);
							// this is my class for defining all relevant kinematic variables

							if (variables.channel_test(variables)) {

								long status = 0;

								int helicity = variables.get_helicity(); 
								// helicity of event, might be 0

								// lab kinematics
								double e_p = variables.e_p();
								double e_theta = variables.e_theta();
								double e_phi = variables.e_phi();
								double p1_p = variables.p1_p();
								double p1_theta = variables.p1_theta();
								double p1_phi = variables.p1_phi();
								double p2_p = variables.p2_p();
								double p2_theta = variables.p2_theta();
								double p2_phi = variables.p2_phi();

								double e_px = variables.e_px();
								double e_py = variables.e_py();
								double e_pz = variables.e_pz();
								double p1_px = variables.p1_px();
								double p1_py = variables.p1_py();
								double p1_pz = variables.p1_pz();
								double p2_px = variables.p2_px();
								double p2_py = variables.p2_py();
								double p2_pz = variables.p2_pz();

								// DIS variables
								double Q2 = variables.Q2();
								double W = variables.W();
								double y = variables.y();
								double Mx = variables.Mx(); // Mx(e:p1:p2:X);
								double Mx1 = variables.Mx1(); // Mx(e:p1:X);
								double Mx2 = variables.Mx2(); // Mx(e:p2:X);
								double Mx3 = variables.Mx3(); // Mx(e:X);

								// SIDIS variables
								double x = variables.x();
								double z = variables.z();
								double xF = variables.xF();
								double pT = variables.pT();
								double eta = variables.eta();
								double eta_gN = variables.eta_gN();

								// SIDIS dihadron variables
								double z1 = variables.z1();
								double z2 = variables.z2();
								double xF1 = variables.xF1();
								double xF2 = variables.xF2();
								double zeta = variables.zeta(); 
								// assumption here is p1 is a proton
								// otherwise zeta is meaningless
								double Mh = variables.Mh();
								double pT1 = variables.pT1();
								double pT2 = variables.pT2();
								double pTpT = variables.pTpT();
								double eta1 = variables.eta1();
								double eta2 = variables.eta2();
								double Delta_eta = variables.Delta_eta();
								double eta1_gN = variables.eta1_gN();
								double eta2_gN = variables.eta2_gN();

								// angles 
								double phi1 = variables.phi1();
								double phi2 = variables.phi2();
								double Delta_phi = variables.Delta_phi();
								double phih = variables.phih();
								double phiR = variables.phiR();
								double theta = variables.theta();

								// vertices 
								double vz_e = variables.vz_e();
								double vz_p1 = variables.vz_p1();
								double vz_p2 = variables.vz_p2();

								// helicity check
								if (helicity == 0 ) { continue; }

								// QA check
								if (!qa.OkForAsymmetry(runnum,evnum)) { status += 1e0; }

								// vertex checks
								if (Math.abs(vz_e + 3) > 5) { status += 1e1; }
								if (Math.abs(vz_e - vz_p1) > 7) { status += 1e2; }
								if (Math.abs(vz_e - vz_p2) > 7) { status += 1e3; }
								if (Math.abs(vz_p1 - vz_p2) > 5) { status += 1e4; }

								// we decided to only dump e_theta > 8 
								// but I processed down to 6 degrees so this continue is necessary
								if (e_theta < 8*3.14159/180) { continue; }
								// but the analysis will be done on e_theta > 10
								if (e_theta < 10*3.14159/180 || 
									e_theta > 30*3.14159/180) { status += 1e5; }
								if (p1_theta > 30*3.14159/180) { status += 1e6; }
								if (p2_theta > 30*3.14159/180) { status += 1e7; }

								// momenta cuts
								if (p1_p < 0.5) { status += 1e8; }
								if (p2_p < 1.2 || p2_p > 4.0) { status += 1e9; }

								// radiative cuts
								if (y > 0.75) { status += 1e10; }

								// z cuts (sort of removing target fragmenation pions)
								if (z2 < 0.2 || z2 > 0.7) { status += 1e11; }

								// current vs target fragmentation cuts
								if (xF1 > 0) { status += 1e12; }
								if (xF2 < 0) { status += 1e13; }
								if (Delta_eta < 0) { status += 1e14; } 

								// missing mass cuts
								if (Mx < 0.85) { status += 1e15; }
								if (Mx1 < 1.6) { status += 1e16; }


								// if (status != 0) { continue; }
								// append event to next line of the text file
								file.append(status+" "+runnum+" "+evnum+" "+helicity+" ");
								file.append(e_p+" "+e_theta+" "+e_phi+" "+vz_e+" ");
								file.append(p2_p+" "+p2_theta+" "+p2_phi+" "+vz_p2+" ");
								file.append(p1_p+" "+p1_theta+" "+p1_phi+" "+vz_p1+" ");
								file.append(Q2+" "+W+" "+x+" "+y+" "+z2+" "+z1+" ");
								file.append(Mx+" "+Mx1+" "+Mx2+" ");
								file.append(pT2+" "+pT1+" "+pTpT+" ");
								file.append(xF2+" "+xF1+" "+eta2+" "+eta1+" "+Delta_eta+" ");
								file.append(phi2+" "+phi1+" "+Delta_phi+"\n");
							}
						}
					}
				}
			}
			println(); println();
			print("1: status, 2: runnum, 3: evnum, 4: helicity, ");
			print("5: e_p, 6: e_theta, 7: e_phi, 8: vz_e, ");
			print("9: pi_p, 10: pi_theta, 11: pi_phi, 12: vz_pi, ");
			print("13: P_p, 14: P_theta, 15: P_phi, 16: vz_P, ");
			print("17: Q2, 18: W, 19: x, 20: y, 21: z_pi, 22: z_P, ");
			print("23: Mx(e:pi:P:X), 24: Mx(e:pi:X), 25: Mx(e:P:X), ");
			print("26: pT_pi, 27: pT_P, 28: pTpT, ");
			print("29: xF_pi, 30: xF_P, 31: eta_pi, 32: eta_P, 33: Delta_eta, ");
			print("34: phi1 (gamma*N COM), 35: phi2 (gamma*N COM), 36: Delta_phi. ");


			println(); println();
			println("Set p1 PID = "+p1_Str+"\n");
			println("Set p2 PID = "+p2_Str+"\n");
			println("output file is: "+file);
			println(); println();
			println("Please be sure your notation of p1 and p2 matches your intention!")
			println("p1 in this code should be a proton -- for things like zeta calculation.")
			println("This is the opposite of the notation of Aram Kotzinian.")
		}

	}
}