package mantesting;

import solver.*;
import ttp.TTP1Instance;
import ttp.TTPSolution;
import utils.Deb;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.concurrent.*;

public class CLI {

  public static void main(String[] args) {

    if (args.length < 2) {
      args = new String[]{"a280_n2790_uncorr_10.ttp"};
    }

    String[] spl = args[0].split("_",2);

    // TTP instance name
    final String inst = args[0];

    // algorithm name
    final String algoName = "savi";

    // output file
    final String outputFile;
    if (args.length >= 2)
      outputFile = "./output/"+args[1];
    else
      outputFile = "./output/"+algoName+".csv";

    // runtime limit
    long runtimeLimit = 600;
    if (args.length >= 3)
      runtimeLimit = Long.parseLong(args[2]);

    // TTP instance
    final TTP1Instance ttp = new TTP1Instance(spl[0]+"-ttp/"+inst);

    /* algorithm to run */
    SearchHeuristic algo = new SAVI(ttp);

    // runnable class
    class TTPRunnable implements Runnable {

      String resultLine;
      TTPSolution sx;

      @Override
      public void run() {
        /* start search & measure runtime */
        long startTime, stopTime;
        long exTime;
        startTime = System.currentTimeMillis();

        sx = algo.search();

        stopTime = System.currentTimeMillis();
        exTime = stopTime - startTime;

        /* print result */
        resultLine = inst + " " + Math.round(sx.ob) + " " + (exTime/1000.0);

      }
    };

    // my TTP runnable
    TTPRunnable ttprun = new TTPRunnable();
    ExecutorService executor = Executors.newFixedThreadPool(4);
    Future<?> future = executor.submit(ttprun);
    executor.shutdown();  // reject all further submissions

    // limit execution time to 600 seconds
    try {
      future.get(runtimeLimit, TimeUnit.SECONDS);  // wait X seconds to finish
    } catch (InterruptedException e) {
      System.out.println("job was interrupted");
    } catch (ExecutionException e) {
      System.out.println("caught exception: " + e.getCause());
    } catch (TimeoutException e) {
      future.cancel(true);
      System.out.println("/!\\ Timeout");
    }

    // wait for execution to be done
    try {
      executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    } catch (InterruptedException e) {
      e.printStackTrace();
    }

    // print results
    Deb.echo(ttprun.resultLine);

    // log results into text file
    try {
      File file = new File(outputFile);
      if (!file.exists()) file.createNewFile();
      Files.write(Paths.get(outputFile), (ttprun.resultLine+"\n").getBytes(), StandardOpenOption.APPEND);
    } catch (IOException e) {
      e.printStackTrace();
    }

    // save solution in a file
    try {
      String currentTime = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss").format(new Date());
      PrintWriter pw = new PrintWriter("./output/solutions/"+inst+"-"+algoName+"-"+currentTime+".txt");
      pw.println(ttprun.sx);
      pw.close();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    }

  }
}
