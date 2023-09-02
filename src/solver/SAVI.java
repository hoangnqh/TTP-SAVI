package solver;

import ttp.TTP1Instance;
import ttp.TTPSolution;
import utils.Deb;
import utils.RandGen;

/**
 Original source code: https://github.com/yafrani/ttplab

 Improved by Hoang Nguyen:
  - function TSPSimulatedAnnealing: Use SA with vertex insertion to optimize tour.
  - function search
 */

public class SAVI extends LocalSearch {

  public double T_abs;       // absolute temperature
  public double T0;           // initial temperature
  public double alpha;       // cooling rate
  public double trialFactor; // number of trials (per temperature)


  public SAVI() {
    super();
    // use default config
    SAConfig();
  }

  public SAVI(TTP1Instance ttp) {
    super(ttp);
    // use default config
    SAConfig();
  }


  // SA params config
  // default config
  void SAConfig() {

    int nbItems = ttp.getNbItems();

    T_abs = 1;
//    T0 = 100.0;
//    alpha = 0.95;
    alpha=0.9578;
    T0=98;

    trialFactor = generateTFLinFit(nbItems);
  }


  /**
   * simulated annealing
   *
   * deal with the KRP sub-problem
   * this function applies a simple bit-flip
   */
  public TTPSolution simulatedAnnealing(TTPSolution sol) {

    // copy initial solution into improved solution
    TTPSolution sBest = sol.clone();

    // TTP data
    int nbCities = ttp.getNbCities();
    int nbItems = ttp.getNbItems();
    int[] A = ttp.getAvailability();
    double maxSpeed = ttp.getMaxSpeed();
    double minSpeed = ttp.getMinSpeed();
    long capacity = ttp.getCapacity();
    double C = (maxSpeed - minSpeed) / capacity;
    double R = ttp.getRent();

    // initial solution data
    int[] tour = sol.getTour();
    int[] pickingPlan = sol.getPickingPlan();

    // delta parameters
    int deltaP, deltaW;

    // best solution
    double GBest = sol.ob;

    // neighbor solution
    long fp;
    double ft, G;
    long wc;
    int origBF;
    int k, r;
    int nbIter = 0;

    double T = T0;
    long trials = Math.round(nbItems*trialFactor);

    if (debug) Deb.echo(">>>> TRIAL FACTOR: "+trialFactor);

    //===============================================
    // start simulated annealing process
    //===============================================
    do {
      nbIter++;

      // cleanup and stop execution if interrupted
      if (Thread.currentThread().isInterrupted()) break;

      for (int u=0; u<trials; u++) {

        // browse items randomly
        k = RandGen.randInt(0, nbItems - 1);

        // check if new weight doesn't exceed knapsack capacity
        if (pickingPlan[k] == 0 && ttp.weightOf(k) > sol.wend) continue;

        // calculate deltaP and deltaW
        if (pickingPlan[k] == 0) {
          deltaP = ttp.profitOf(k);
          deltaW = ttp.weightOf(k);
        } else {
          deltaP = -ttp.profitOf(k);
          deltaW = -ttp.weightOf(k);
        }
        fp = sol.fp + deltaP;

        // handle velocity constraint
        // index where Bit-Flip happened
        origBF = sol.mapCI[A[k] - 1];
        // starting time
        ft = origBF == 0 ? .0 : sol.timeAcc[origBF - 1];
        // recalculate velocities from bit-flip city
        // to recover objective value
        for (r = origBF; r < nbCities; r++) {
          wc = sol.weightAcc[r] + deltaW;
          ft += ttp.distFor(tour[r] - 1, tour[(r + 1) % nbCities] - 1) / (maxSpeed - wc * C);
        }
        // compute recovered objective value
        G = fp - ft * R;

        //=====================================
        // update if improvement or
        // Boltzmann condition satisfied
        //=====================================
        double mu = Math.random();
        double energy_gap = G - GBest;
        boolean acceptance = energy_gap > 0 || Math.exp(energy_gap / T) > mu;
        if (acceptance) {

          GBest = G;

          // bit-flip
          pickingPlan[k] = pickingPlan[k] != 0 ? 0 : A[k];

          //===========================================================
          // recover accumulation vectors
          //===========================================================
          if (pickingPlan[k] != 0) {
            deltaP = ttp.profitOf(k);
            deltaW = ttp.weightOf(k);
          } else {
            deltaP = -ttp.profitOf(k);
            deltaW = -ttp.weightOf(k);
          }
          fp = sol.fp + deltaP;
          origBF = sol.mapCI[A[k] - 1];
          ft = origBF == 0 ? 0 : sol.timeAcc[origBF - 1];
          for (r = origBF; r < nbCities; r++) {
            // recalculate velocities from bit-flip city
            wc = sol.weightAcc[r] + deltaW;
            ft += ttp.distFor(tour[r] - 1, tour[(r + 1) % nbCities] - 1) / (maxSpeed - wc * C);
            // recover wacc and tacc
            sol.weightAcc[r] = wc;
            sol.timeAcc[r] = ft;
          }
          G = fp - ft * R;
          sol.ob = G;
          sol.fp = fp;
          sol.ft = ft;
          sol.wend = capacity - sol.weightAcc[nbCities - 1];
          //===========================================================

        }

      }

      // update best if improvement
      if (sol.ob > sBest.ob) {
        sBest = sol.clone();
      }

      if (this.debug) {
        Deb.echo(">> KRP " + nbIter + ": ob=" +
                String.format("%.0f",sol.ob));
      }

      // cool down temperature
      T = T * alpha;

      // stop when temperature reach absolute value
    } while (T > T_abs);


    // in order to recover all history vector
    ttp.objective(sBest);

    return sBest;
  }

  public TTPSolution TSPSimulatedAnnealing(TTPSolution sol) {
    // copy initial solution into improved solution
    TTPSolution sBest = sol.clone();
    long profitFinal = sol.fp;

    // TTP data
    int nbCities = ttp.getNbCities();
    int nbItems = ttp.getNbItems();
    int[] A = ttp.getAvailability();
    double maxSpeed = ttp.getMaxSpeed();
    double minSpeed = ttp.getMinSpeed();
    long capacity = ttp.getCapacity();
    double C = (maxSpeed - minSpeed) / capacity;
    double R = ttp.getRent();

    // initial solution data
    int[] tour = sol.getTour();
    int[] pickingPlan = sol.getPickingPlan();

    long [] weightRec = new long[nbCities];
    for (int i = 0; i < nbCities; i++){
      int city = tour[i] - 1;
      weightRec[city] = sol.weightRec[i];
    }

    // delta parameters
    int deltaP, deltaW;

    // best solution
    double GBest = sol.ob;

    int c1, c2, c3, c4;
    long weight, c1Weight, c2Weight, c3Weight;
    int nbIter = 0;

    double T = T0;
    // 6762000
    int trials = Math.min(nbCities * 10, (12762000 / nbCities));
//    int trials = nbCities;
//    System.out.println("nbCites: "+nbCities+", trials: "+trials);




    //===============================================
    // start simulated annealing process
    //===============================================
    do {
      nbIter++;

      for (int u=0; u<trials; u++) {
        // cleanup and stop execution if interrupted
        if (Thread.currentThread().isInterrupted()) break;

        int pos_i = RandGen.randInt(1, nbCities - 1);

        int[] newTour = sol.getTour().clone();

        int pos = pos_i;
        while (pos > 1) {
          int temp = newTour[pos];
          newTour[pos] = newTour[pos - 1];
          newTour[pos - 1] = temp;
          pos -= 1;
        }

        TTPSolution newSolution = new TTPSolution(newTour, pickingPlan);
//        ttp.objective(newSolution);
        long new_weight = 0;
        newSolution.ft = 0;
        for (int i = 0; i < nbCities; i++){
          int city = newTour[i] - 1;
          int nextCity = newTour[(i + 1) % nbCities] - 1;
          new_weight += weightRec[city];
//          if (city == 0) System.out.println("weightRec = "+weightRec[city]);
          newSolution.weightRec[i] = weightRec[city];
          newSolution.weightAcc[i] += (i == 0 ? 0 : newSolution.weightAcc[i - 1]) + weightRec[city];
          newSolution.ft += ttp.distFor(city, nextCity) / (maxSpeed - new_weight * C);
//          newSolution.weightAcc[i] = new_weight;
        }
        newSolution.ob = profitFinal - R * newSolution.ft;
        double tempBestZ = newSolution.ob;
        double posBestZ = 1;
        double[] time = new double[nbCities];
//        int[] newTour = newTour.clone();

        time[1] = newSolution.ft;
        weight = 0;
        for (int j = 2; j < nbCities; j++){
          time[j] = time[j - 1];

          c1 = (j == 2 ? 0 : j - 1);
          c1 = j - 2;

          c2 = 1;
          c2 = j - 1;

          c3 = j;
          c3 = j;

          c4 = (j + 1) % nbCities;
//          System.out.println(j);
//          System.out.println(c1+" "+c2+" "+c3+" "+c4);

          // UPDATE
          // Lost time
          // before: c1 c2 c3 c4
          // c1 -> c2
//          c1Weight = (c1 > 0 ? (newSolution.weightAcc[c1] - newSolution.weightAcc[c1 - 1] + weight) : 0);
          c1Weight =  weightRec[newTour[c1] - 1] + weight;
          time[j] -= ttp.distFor(newTour[c1] - 1, newTour[c2] - 1) / (maxSpeed - c1Weight * C);
          // c2 -> c3
//          c2Weight = newSolution.weightAcc[c2] - newSolution.weightAcc[c2 - 1] + c1Weight;
          c2Weight =  weightRec[newTour[c2] - 1] + c1Weight;
          time[j] -= ttp.distFor(newTour[c2] - 1, newTour[c3] - 1) / (maxSpeed - c2Weight * C);
          // c3 -> c4
//          c3Weight = newSolution.weightAcc[c3] - newSolution.weightAcc[c3 - 1] + c2Weight;
          c3Weight =  weightRec[newTour[c3] - 1] + c2Weight;
          time[j] -= ttp.distFor(newTour[c3] - 1, newTour[c4] - 1) / (maxSpeed - c3Weight * C);


          // Extra time
          // after:  c1 c3 c2 c4
          // c1 -> c3
          time[j] += ttp.distFor(newTour[c1] - 1, newTour[c3] - 1) / (maxSpeed - c1Weight * C);
          // c3 -> c2
//          c3Weight = newSolution.weightAcc[c3] - newSolution.weightAcc[c3 - 1] + c1Weight;
          c3Weight = weightRec[newTour[c3] - 1] + c1Weight;
          time[j] += ttp.distFor(newTour[c3] - 1, newTour[c2] - 1) / (maxSpeed - c3Weight * C);
          // c2 -> c4
//          c2Weight = newSolution.weightAcc[c2] - newSolution.weightAcc[c2 - 1] + c3Weight;
          c2Weight = weightRec[newTour[c2] - 1] + c3Weight;
          time[j] += ttp.distFor(newTour[c2] - 1, newTour[c4] - 1) / (maxSpeed - c2Weight * C);

          // update c1Weight after swap(c2, c3)
          weight = c1Weight;

          int temp_ = newTour[j - 1];
          newTour[j - 1] = newTour[j];
          newTour[j] = temp_;

          double presentZ = profitFinal - ttp.getRent() * time[j];
          if (tempBestZ < presentZ){
            tempBestZ = presentZ;
            posBestZ = j;
          }
        }

        double mu = Math.random();
        double energy_gap = tempBestZ - GBest;
        boolean acceptance = energy_gap > 0 || Math.exp(energy_gap / T) > mu;
        if (acceptance) {
//          newTour = sBest.getTour().clone();
          newTour = sol.getTour().clone();
          pos = pos_i;
          while (pos > posBestZ) {
            int temp = newTour[pos];
            newTour[pos] = newTour[pos - 1];
            newTour[pos - 1] = temp;
            pos -= 1;
          }
          while (pos < posBestZ) {
            int temp = newTour[pos];
            newTour[pos] = newTour[pos + 1];
            newTour[pos + 1] = temp;
            pos += 1;
          }
          GBest = tempBestZ;
          sol = new TTPSolution(newTour, pickingPlan);
//           ttp.objective(sol);
          new_weight = 0;
          sol.ft = 0;
          for (int i = 0; i < nbCities; i++){
            int city = newTour[i] - 1;
            int nextCity = newTour[(i + 1) % nbCities] - 1;
            new_weight += weightRec[city];
            sol.ft += ttp.distFor(city, nextCity) / (maxSpeed - new_weight * C);
            sol.weightAcc[i] = new_weight;
          }
          sol.ob = profitFinal - R * sol.ft;
        }
      }

      // update best if improvement
      if (sol.ob > sBest.ob) {
//        System.out.println("TSPSA: "+sol.ob);
        sBest = sol.clone();
        //solution = ttp.evaluate(tour, pickingPlan, true);
      }

      // cool down temperature
      T = T * alpha;

      // stop when temperature reach absolute value
    } while (T > T_abs);

    // in order to recover all history vector
    ttp.objective(sBest);

    return sBest;
  }

  @Override
  public TTPSolution search() {
    //===============================================
    // generate initial solution
    //===============================================
    if (s0==null) {
      Constructive construct = new Constructive(ttp);
      // use Lin-Kernighan to initialize the tour
      s0 = new TTPSolution(
              construct.linkernTour(),
              construct.zerosPickingPlan()
      );

      // pre-process the knapsack
      // insert and eliminate items
      if (ttp.getNbCities() < 30000) s0 = insertAndEliminate(s0);
      else s0 = insertT2(s0);
      ttp.objective(s0);
      //System.out.println(s0.ob);
//      s0 = lsBitFlip(s0);

//      Initialization init = new Initialization(ttp);
//      s0 = init.lkPackIterative();
    }
    ttp.objective(s0);
    if (this.debug) {
      Deb.echo("STARTING SOL >> " + s0.ob);
    }
    //===============================================

    // copy initial solution into improved solution
    TTPSolution sol = s0.clone();
//    System.out.println("Start: "+String.format("%.2f", sol.ob));

    // best found
    double GBest = sol.ob;
    // number of iterations
    int nbIter = 0;
    // improvement tag
    boolean improved;
    //===============================================
    // start cosolver search
    //===============================================

    double mark1 = -Double.MAX_VALUE;
    double mark2 = mark1;
    double mark3 = mark1;

    do {
      nbIter++;
      improved = false;

      // 2-opt heuristic on TSKP
      if (mark1 == sol.ob) break;
      sol = fast2opt(sol);
      mark1 = sol.ob;
      //System.out.println("2opt:  "+String.format("%.2f", sol.ob));

      // stop execution if interrupted
      if (Thread.currentThread().isInterrupted()) return sol;
      if (mark2 == sol.ob) break;
      sol = TSPSimulatedAnnealing(sol);
      mark2 = sol.ob;
      //System.out.println("TSPSA: "+String.format("%.2f", sol.ob));

      // stop execution if interrupted
      if (Thread.currentThread().isInterrupted()) return sol;
      // simple bit-flip on KRP
      if (mark3 == sol.ob) break;
      sol = simulatedAnnealing(sol);
      mark3 = sol.ob;
      //System.out.println("SA:    "+String.format("%.2f", sol.ob));

      // update best if improvement
      if (sol.ob > GBest) {
        GBest = sol.ob;
        improved = true;
      }
      // stop execution if interrupted
      if (Thread.currentThread().isInterrupted()) return sol;

      // debug msg
      if (this.debug) {
        Deb.echo("Best "+nbIter+":");
        Deb.echo("ob-best: "+sol.ob);
        Deb.echo("wend   : "+sol.wend);
        Deb.echo("---");
      }

      // stop when no improvements
    } while (improved);
    //===============================================
    sol = lsBitFlip(sol);
    return sol;
  }



  public static double generateTFLinFitx(int xi) {
    int i=-1;
//    int[] x = new int[]{         50,    204,  609,  1147,  8034, 38875, 105318,  253568, 338090};
    int[] x = new int[]{         1,   75,    375,  790,  2102, 15111, 70250, 140500,  338090};
    double[] y = new double[]{57872,  13896,  700,   350,    16,     1,   0.16,  0.0493,   0.03};

//    int[] x = new int[]{         1,   75,    375,  790,  2102, 15111, 70250, 140500,  338090};
//    double[] y = new double[]{
//      57872, 13896, 350, 16, 1, .16, .0493, .03
//    };

    int n=y.length;
    for (int k=0; k<n-1; k++) {
      if (x[k] <= xi && xi < x[k + 1]) {
        i = k;
        break;
      }
    }
    if (xi <= x[0]) {
      return 57872.0;
    }
    if (xi >= x[n-1]) {
      return 0.03;
    }

    double m = ( y[i]-y[i+1] ) / ( x[i]-x[i+1] );
    double b = y[i]-m*x[i];

    double yi = m*xi + b;

    return yi;
  }

  public static double generateTFLinFit(int xi) {
    int i=-1;
//    int[] z = new int[]{ 1,  75, 375, 790, 2102, 15111, 70250, 140500, 338090};

//    int[] x = new int[]{         1,   75,    375,  790,  2102, 15111, 70250, 140500,  338090};
//    double[] y = new double[]{57872,  13896,  700,   350,    16,     1,   0.16,  0.0493,   0.03};

    int[] x = new int[]{ 1, 130, 496, 991, 3038, 18512, 75556, 169046, 338090};
    double[] y = new double[]{57872,  13896,  700,   350,    16,     1,   0.16,  0.0493,   0.03};


    int n=y.length;
    for (int k=0; k<n-1; k++) {
      if (x[k] <= xi && xi < x[k + 1]) {
        i = k;
        break;
      }
    }
    if (xi <= x[0]) {
      return 57872.0;
    }
    if (xi >= x[n-1]) {
      return 0.03;
    }

    double m = ( y[i]-y[i+1] ) / ( x[i]-x[i+1] );
    double b = y[i]-m*x[i];

    double yi = m*xi + b;

    return yi;
  }

  public static double generateTFExpFit(int X) {
    double a =   1.614e+05;
    double b =   -0.006915;
    double c =       145.4;
    double d =  -5.191e-05;

    double A = a * Math.exp(b * X) + c * Math.exp(d * X);

    return A;
  }

  public static double generateTFManualFit(int X) {
    int nbItems = X;
    double trialFactor;
    if (nbItems < 500)
      trialFactor = 1000;// reduce... (time)

    else if (nbItems < 1000)
      trialFactor = 100;
    else if (nbItems < 5000)
      trialFactor = 50;

    else if (nbItems < 20000)
      trialFactor = 10; //was 5... retest others
    else if (nbItems < 100000)
      trialFactor = 1;
    else if (nbItems < 200000)
      trialFactor = .04;
    else
      trialFactor = .03;
    return trialFactor;
  }

}