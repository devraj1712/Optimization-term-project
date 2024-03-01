#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

double conjugate_grad(double x[],double uplim[], double lowlim[],int n,int m, double R,int fevalcg);
double bounding_phase(double x[], double s[], double limits[], int n, double R);
double newton_raphson(double xprev, double x0,double x[],double s[],int feval2, int n, double R);
double obj_func(double x[]);

int main()
{
    int i,m,n,num,t,c,feval_seq,feval_total;
    FILE *file1=fopen("q1_input.txt","r");
    fscanf(file1,"no_of_variables %d\n",&n);
    fscanf(file1,"no_of_constraints %d\n",&m);
    int uint[n],lint[n];
    double x[n],lowlim[n],uplim[n],R,epsilon,P,P_prev,g[m],f;
    for(i=0;i<n;i++)
    {
        fscanf(file1,"lower_limit %*s %lf\n",&lowlim[i]);
        fscanf(file1,"upper_limit %*s %lf\n",&uplim[i]);
        printf("%lf \t %lf \n",lowlim[i],uplim[i]);
    }
    fclose(file1);
    srand(time(0));
    for(i=0;i<n;i++)
    {
        uint[i]=floor(uplim[i]);  // Randomly generating initial guess value between limits read from file
        lint[i]=floor(lowlim[i]);
        num=rand()%(uint[i]-lint[i]);
        x[i]=num+lowlim[i];
        printf("%lf \t",x[i]);
    }
    t=0;
    c=10;
    R=2;
    epsilon=0.001;
    feval_total=0;
    do
    {
        if(t>0)
        {
            P_prev=P;
        }
        else
        {
            P_prev=0;
        }
        printf("\n Initializing conjugate gradient method taking R equal to %lf \n",R);
        printf("sequence no. %d\n",t);
        P=conjugate_grad(x,uplim,lowlim,n,m,R,feval_seq);
        R=c*R;
        t=t+1;
        feval_total=feval_total+feval_seq;
    }while(sqrt( pow(P-P_prev,2)>epsilon) && R<pow(10,6) );
    f=obj_func(x);
    printf("\nAfter sequence number %d f(x) = %lf \n",t,f);
    return 0;
}

double pen_fun(double x[], int n, double R,double g[])
{
    //NOTE - Change input file name when changing objective function
    FILE *of = fopen("constraint violations.txt","a");

    //1. Problem 1
    fprintf(of,"R \t g1 \t g2 \n");
    double gnorm[2],g1max,g2max;
    g[0] = pow(x[0]-5,2) + pow(x[1]-5,2) - 100;
    g[1] = 82.81 - pow(x[0]-6,2) - pow(x[1]-5,2);
    if(g[0]>=0) // bracket operator returns only negative values
    {
        g[0]=0;
    }
    if(g[1]>=0)
    {
        g[1]=0;
    }
    //normalizing the constraints
    g1max = -35;
    g2max = -138.19;
    //printf("\nConstraint violation g1 = % lf \n",g1);
    //printf("Constraint violation g2 = % lf \n",g2);
    fprintf(of,"%lf \t %lf \t %lf \n",R,g[0],g[1]);
    fclose(of);
    gnorm[0]=g[0]/g1max;
    gnorm[1]=g[1]/g2max;
    return pow(x[0]-10,3) + pow(x[1]-20,3) + R*(pow(gnorm[0],2) + pow(gnorm[1],2) );

    //2. Problem 2
    /*
    fprintf(of,"R \t g1 \t g2 \n");
    double gnorm[2],g1max,g2max;
    g[0] = -1*pow(x[0],2) + x[1] - 1;
    g[1] = -1 + x[0] - pow(x[1]-4,2);
    if(g[0]>=0) // bracket operator returns only negative values
    {
        g[0]=0;
    }
    if(g[1]>=0)
    {
        g[1]=0;
    }
    //normalizing the constraints
    g1max = -101;
    g2max = -37;
    // printf("\nConstraint violation g1 = % lf \n",g1);
    // printf("Constraint violation g2 = % lf \n",g2);
    fprintf(of,"%lf \t %lf \t %lf \n",R,g[0],g[1]);
    fclose(of);
    gnorm[0]=g[0]/g1max;
    gnorm[1]=g[1]/g2max;
    return -1*(pow( sin(2*3.142857*x[0]),3 )*sin(2*3.142857*x[1]))/( pow(x[0],3)*(x[0] + x[1]) ) + R*(pow(gnorm[0],2) + pow(gnorm[1],2) );
    */
    //3. Problem 3
    /*
    fprintf(of,"R \t g1 \t g2 \t g3 \t g4 \t g5 \t g6\n");
    double g1max,g2max,g3max,g4max,g5max,g6max,gnorm[6];
    g[0] = 1 - 0.0025*(x[3] + x[5]);
    g[1] = 1 - 0.0025*(-1*x[3] + x[4] + x[6]);
    g[2] = 1 - 0.01*(-1*x[5] + x[7]);
    g[3] = -100*x[0] + x[0]*x[5] - 833.33252*x[3] + 83333.333;
    g[4] = x[1]*(-1*x[3] + x[6]) + 1250*(x[3] - x[4]);
    g[5] = -1*x[2]*x[4] + x[2]*x[7] + 2500*x[4] - 1250000;
    for(int j=0;j<6;j++)
    {
        if(g[j]>=0) // bracket operator returns only negative values
        {
            g[j]=0;
        }
    }
    //normalizing the constraints
    g1max = -4;
    g2max = -3.975;
    g3max = -8.9;
    g4max = -1748999.187;
    g5max = -11137500;
    g6max = -11125000;
    fprintf(of,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t *lf \n",R,g[0],g[1],g[2],g[3],g[4],g[5]);
    fclose(of);
    gnorm[0]=g[0]/g1max;
    gnorm[1]=g[1]/g2max;
    gnorm[2]=g[2]/g3max;
    gnorm[3]=g[3]/g4max;
    gnorm[4]=g[4]/g5max;
    gnorm[5]=g[5]/g6max;
    return x[0] + x[1] + x[2] + R*(pow(gnorm[0],2) + pow(gnorm[1],2) + pow(gnorm[2],2) + pow(gnorm[3],2) + pow(gnorm[4],2) + pow(gnorm[5],2) );
    */
}

double single_var_pen_func(double a, double x[],double s[],int n, double R)
{
    //NOTE - Change input file name when changing objective function

    //1. Problem 1
    double g1,g2,g1max,g2max;
    g1 = pow((x[0]+a*s[0])-5,2) + pow((x[1]+a*s[1])-5,2) - 100;
    g2 = 82.81 - pow((x[0]+a*s[0])-6,2) - pow((x[1]+a*s[1])-5,2);
    if(g1>=0)  // bracket operator returns only negative values
    {
        g1=0;
    }
    if(g2>=0)
    {
        g2=0;
    }
    g1max = -35;
    g2max = -138.19;
    g1=g1/g1max;
    g2=g2/g2max;

    return pow((x[0]+a*s[0])-10,3) + pow((x[1]+a*s[1])-20,3) + R*( pow(g1,2) + pow(g2,2) );

    //2. Problem 2
    /*
    double g1,g2,g1max,g2max;
    g1 = -1*pow((x[0]+a*s[0]),2) + (x[1]+a*s[1]) - 1;
    g2 = -1 + (x[0]+a*s[0]) - pow((x[1]+a*s[1])-4,2);
    if(g1>=0) // bracket operator returns only negative values
    {
        g1=0;
    }
    if(g2>=0)
    {
        g2=0;
    }
    //normalizing the constraints
    g1max = -101;
    g2max = -37;
    g1=g1/g1max;
    g2=g2/g2max;
    return -1*(pow( sin(2*3.142857*(x[0]+a*s[0])),3 )*sin(2*3.142857*(x[1]+a*s[1])))/( pow((x[0]+a*s[0]),3)*((x[0]+a*s[0]) + (x[1]+a*s[1])) ) + R*(pow(g1,2) + pow(g2,2) );
    */
    //3. Problem 3
    /*
    double g[6],g1max,g2max,g3max,g4max,g5max,g6max;
    g[0] = 1 - 0.0025*((x[3]+a*s[3]) + (x[5]+a*s[5]));
    g[1] = 1 - 0.0025*(-1*(x[3]+a*s[3]) + (x[4]+a*s[4]) + (x[6]+a*s[6]));
    g[2] = 1 - 0.01*(-1*(x[5]+a*s[5]) + (x[7]+a*s[7]));
    g[3] = -100*(x[0]+a*s[0]) + (x[0]+a*s[0])*(x[5]+a*s[5]) - 833.33252*(x[3]+a*s[3]) + 83333.333;
    g[4] = (x[1]+a*s[1])*(-1*(x[3]+a*s[3]) + (x[6]+a*s[6])) + 1250*((x[3]+a*s[3]) - (x[4]+a*s[4]));
    g[5] = -1*(x[2]+a*s[2])*(x[4]+a*s[4]) + (x[2]+a*s[2])*(x[7]+a*s[7]) + 2500*(x[4]+a*s[4]) - 1250000;
    for(int j=0;j<6;j++)
    {
        if(g[j]>=0) // bracket operator returns only negative values
        {
            g[j]=0;
        }
    }
    //normalizing the constraints
    g1max = -4;
    g2max = -3.975;
    g3max = -8.9;
    g4max = -1748999.187;
    g5max = -11137500;
    g6max = -11125000;
    g[0]=g[0]/g1max;
    g[1]=g[1]/g2max;
    g[2]=g[2]/g3max;
    g[3]=g[3]/g4max;
    g[4]=g[4]/g5max;
    g[5]=g[5]/g6max;
    return (x[0]+a*s[0]) + (x[1]+a*s[1]) + (x[2]+a*s[2]) + R*(pow(g[0],2) + pow(g[1],2) + pow(g[2],2) + pow(g[3],2) + pow(g[4],2) + pow(g[5],2) );
    */
}

double obj_func(double x[])
{
    //1. Problem 1
    return pow(x[0]-10,3) + pow(x[1]-20,3) ;

    //2. Problem 2
    //return -1*(pow( sin(2*3.142857*x[0]),3 )*sin(2*3.142857*x[1]))/( pow(x[0],3)*(x[0] + x[1]) ) ;

    //3. Problem 3
    //return x[0] + x[1] + x[2] ;
}
double bounding_phase(double x[], double s[], double limits[], int n, double R)
{
    int i,num,uint,lint;
    double x0,x1,x2,x3,f1,f2,f3,a,b;
    if(limits[0]<limits[1])
    {
        a=limits[0];
        b=limits[1];
    }
    else
    {
        b=limits[0];
        a=limits[1];
    }
    FILE *out_bp;
    printf("\n---------------------------------------------------------------------------------\n");
    printf("Bounding Phase Method\n");
    //Step 1
    uint=floor(b);
    lint=floor(a);
    srand(time(0));
    if(uint==lint) // If the limits for generating random no. are both truncated to same integer, adding an exception
    {
        x0 = (a+b)/2;
    }
    else
    {
        num=rand()%(uint-lint);
        if(num == 0)   // If the random no. generated is a boundary point, adding an exception
        {
            x0 = (a+b)/2;
        }
        else
        {
            x0 = num + a;
        }
    }
    double delta = (b-a)/10;
    int fevalbp;// Function Evaluations
    int k=0;// No. of Iterations
    //Step 2
    printf("\n a %lf \t b %lf \t x0 %lf \n",a,b,x0);
    while(1)
    {   fevalbp=0;
        x1=x0-delta; //New Points
        x2=x0;
        x3=x0+delta;
        f1=single_var_pen_func(x1,x,s,n,R); // Calculate Objective Function
        f2=single_var_pen_func(x2,x,s,n,R);
        f3=single_var_pen_func(x3,x,s,n,R);
        printf("\n x1 %lf \t x3 %lf \t f1 %lf \t f3 %lf \n",x1,x3,f1,f3);
        fevalbp=fevalbp+3;
        if(f1>=f2 && f2>=f3) // Check whether delta should be taken +ve or -ve
        {
            printf("Delta is positive");
            break;
        }
        else if(f1<=f2 && f2<=f3)
        {
            printf("Delta is negative");
            delta=delta*(-1);
            break;
        }
        else
        {
            printf("\n since f2<f1 & f2<f3 , choosing new initial guess for alpha\n");
            x0=(x0+a)/2;
            delta = (b-a)/50;
            printf("\n x0 %lf \n",x0);
        }

    }
    out_bp= fopen("Bounding_Phase_Iterations.out","w"); //Output File
    fprintf(out_bp,"#It\t x\t F(x)\n");
    fprintf(out_bp,"%d\t %lf\t %lf\n", k,x,f2);

    //Step 3
    double fold,fnew,xprev;
    do
    {   if(k==0)
        {
            fold=f2; // Storing previous iteration values to print range
            xprev=x0;
        }
        else
        {
            fold=fnew;
            xprev=x0-((pow(2.0,k-1))*delta);
        }
        x0=x0+(pow(2.0,k)*delta); // Update Value of x
        k=k+1;
        //STEP 4
        printf("\n x0 %lf xprev %lf \n",x0,xprev);
        if(xprev < a)
        {
            xprev=a;
        }
        if(x0 < a)
        {
            x0=a;
            break;
        }
        if(xprev > b)
        {
            xprev=b;
        }
        if(x0 > b)
        {
            x0=b;
            break;
        }
        fnew=single_var_pen_func(x0,x,s,n,R);
        printf("\n k %d\t x0 %lf \t %lf \n", k,x0,fnew);
        fevalbp++;
        fprintf(out_bp,"%d\t %lf\t %lf\n", k,x0,fnew);
    }while(fnew<fold);//Termination Condition
    printf("\n---------------------------------------------------------------------------------\n");
    printf("The minimum lies in (%lf , %lf)", xprev,x0);
    printf("\n#Total number of function evaluations: %d",fevalbp);
    // Store in the file
    fprintf(out_bp,"The minimum lies in (%lf , %lf)", xprev,x0);
    fprintf(out_bp,"\n#Total number of function evaluations: %d",fevalbp);
    fclose(out_bp);
    limits[0]=xprev;
    limits[1]=x0;
    return fevalbp;
}

double newton_raphson(double xprev, double x0,double x[],double s[],int feval2,int n,double R)
{
    int uint,lint,ran_num;
    double y,l_lim,u_lim,error,chk,del,y11,y12,f10,f11,f12,first_der,sec_der;
    FILE *out_nr;
    printf("\n---------------------------------------------------------------------------------\n");
    printf("Newton Raphson Method\n");
    if(xprev<x0)
    {
        l_lim=xprev;
        u_lim=x0;
    }
    else
    {
        u_lim=xprev;
        l_lim=x0;
    }
    //Step 1
    uint=floor(x0);
    lint=floor(xprev);
    srand(time(0));

    if(uint==lint) // If the limits for generating random no. are both truncated to same integer, adding an exception
    {
        y = (x0+xprev)/2;

    }
    else
    {
        ran_num=rand()%(uint-lint);
        //printf("here\n");
        if(ran_num == 0)  // If the random no. generated is a boundary point, adding an exception
        {
            y = (x0+xprev)/2;
        }
        else
        {
            y = ran_num + l_lim;
        }
    }
    printf("\n x0 %lf \t xprev %lf \t y %lf \n",x0,xprev,y);
    error=0.001;
    int k=1; // Number of iterations
    int fevalnr=0; // Function Evaluations
    del=(u_lim-l_lim)/100; // Delta value for central difference
    out_nr= fopen("Newton_Raphson_Iterations.out","w"); // Output file
    fprintf(out_nr,"#It\t x\t F\'(x)\t F\"(x)\n");

    do
    {
        y11=y+del; //Points for calculating central difference
        y12=y-del;
        f10=single_var_pen_func(y,x,s,n,R); // Function Values at new points
        f11=single_var_pen_func(y11,x,s,n,R);
        f12=single_var_pen_func(y12,x,s,n,R);
        first_der=(f11-f12)/(2*del); // Calculate First Derivative
        //Step 2
        sec_der=(f11+f12-(2*f10))/pow(del,2.0); // Calculate Second Derivative
        //Step 3
        y=y-(first_der/sec_der); // Update value of y
        //Step 4
        if(first_der<0)  // Modulus of first derivative to check termination condition
        {
            chk=first_der*(-1);
        }
        else
        {
            chk=first_der;
        }
        //Boundary violations
        if(y<l_lim)
        {
            printf("\n y %lf \n",y);
            if(f10>=0)
            {
                y=l_lim;
            }
            else
            {
                y=u_lim;
            }
            break;
        }
        else if (y>u_lim)
        {
            printf("\n y %lf \n",y);
            if(f10>=0)
            {
                y=u_lim;
            }
            else
            {
                y=l_lim;
            }
            break;
        }
        fprintf(out_nr,"%d\t %lf\t %lf\t %lf\n", k,y,first_der, sec_der);
        //printf("\n y %lf \t fdash %lf \n",y,chk);
        k=k+1;
        fevalnr=fevalnr+3;
    }while(chk>error && k<100);//Termination Condition
    printf("\n---------------------------------------------------------------------------------\n");
    if(y == l_lim || y == u_lim)
    {
        printf("\n Boundary violation - Optimum value lies at the boundary point. \n");
    }
    else
    {
        printf("After %d iterations f\'a = %lf < %lf\n", k,first_der,error);
    }
    printf("The optimum value of a = %lf\n",y);
    printf("\n#Total number of function evaluations: %d",fevalnr);
    // Store in the file
    fprintf(out_nr,"The optimum value of a = %lf\n", y);
    fprintf(out_nr,"\n#Total number of function evaluations: %d",fevalnr);
    fclose(out_nr);
    feval2=feval2+fevalnr;
    return y;
}

double conjugate_grad(double x[],double uplim[],double lowlim[],int n,int m, double R, int fevalcg)
{
    int i,iter;
    double error_1,error_2,error_3;
    FILE *out_cg;
    double x1[n],x2[n],grad_fx[n],f0,f1[n],f2[n],s[n],unitS_old[n],unitS_new[n],x_old[n],g[m];
    double del,temp_var,normgrad_sqold,normgrad_sqnew,s_norm,normx_diffsq,normx_sq,dot_prod,angle,alpha;
    int feval2,restart_pt;
    //Step 1
    error_1=error_2=error_3=0.001;
    iter=0;
    restart_pt=0;
    fevalcg=0; // Function Evaluations
    del=0.001; // Delta value for central difference
    printf("\n---------------------------------------------------------------------------------\n");
    printf("Conjugate Gradient Method\n");
    printf("R Value for present sequence is %lf, x0 is : ",R);
    out_cg= fopen("Conjugate_Gradient_Iterations.txt","a"); //Output File
    fprintf(out_cg,"R value for present sequence is %lf, x0 is : ",R);
    for(int i=0;i<n;i++)
    {
        printf("%lf \t",x[i]);
        fprintf(out_cg,"%lf \t",x[i]);
    }
    f0=pen_fun(x,n,R,g); // Function Value at new points
    normgrad_sqnew=0;   //initializing variable to store norm of gradient squared
    fprintf(out_cg,"\n#It\t F(x)\t ");
    for(i=0;i<n;i++)
    {
        fprintf(out_cg,"x%d \t",i+1);
        fprintf(out_cg,"F\'(x%d) \t",i+1);
    }
    for(i=0;i<m;i++)
    {
        fprintf(out_cg,"g%d \t",i+1);
    }
    fprintf(out_cg,"Angle\t||grad_f(x)||\t||x(k)-x(k-1)||/||(x(k-1)||\n");
    do
    {   //Step 2
        for(i=0;i<n;i++)
        {
            x1[i]=x[i]+del; //Points for calculating central difference
            x2[i]=x[i]-del;
        }
        for(i=0;i<n;i++)
        {
            temp_var=x[i];
            x[i]=x1[i];        // Replacing x(i) with x(i)+del keeping other n-1 variables unchanged
            f1[i]=pen_fun(x,n,R,g);   // Computing f(x(i)+del)
            x[i]=x2[i];        // Replacing x(i) with x(i)-del keeping other n-1 variables unchanged
            f2[i]=pen_fun(x,n,R,g);   // Computing f(x(i)-del)
            x[i]=temp_var;
        }
        fevalcg=fevalcg+(2*n);
        normgrad_sqold=normgrad_sqnew;
        normgrad_sqnew=0;
        s_norm=0;
        printf("\n x(%d) \t grad_f(x(%d)) \t s(%d) : \n",iter,iter,iter);
        for(i=0;i<n;i++)
        {
            grad_fx[i]=(f1[i]-f2[i])/(2*del); // Calculate gradient
            normgrad_sqnew =normgrad_sqnew + (grad_fx[i]*grad_fx[i]);  // Compute norm of gradient squared
        }
        for(i=0;i<n;i++)    // Calculate search direction
        {
            if(iter <1)
            {
                s[i]=-1*grad_fx[i];
                s_norm= s_norm + (s[i]*s[i]);   //Compute norm of s0 for checking linear independency in next iteration
            }
            else
            {   //Step 4

                s[i]=-1*grad_fx[i]+ (normgrad_sqnew/normgrad_sqold)*s[i]; //Compute search direction s(k)
                //printf("\n normnew %lf \t normold %lf \t s %lf \n",normgrad_sqnew,normgrad_sqold,s[i]);
                s_norm= s_norm + (s[i]*s[i]);
            }
            printf("%lf \t %lf \t %lf \n",x[i],grad_fx[i],s[i]);
        }
        s_norm=sqrt(s_norm);
        for(i=0;i<n;i++)
        {
            unitS_old[i]=unitS_new[i];
            unitS_new[i]=s[i]/s_norm;
        }
        dot_prod=0;
        if(iter>0)
        {
            for(i=0;i<n;i++)
            {
                dot_prod=dot_prod + (unitS_old[i]*unitS_new[i]);
            }
            angle=acos(dot_prod);
            angle=angle*(180/3.142857); // Computing value of angle between unit direction vectors in degrees
        }
        if(iter>0 && angle<5)  // Check for linear independency
        {
            printf("\nAt iteration no. %d , angle between unit direction vectors is %lf degrees\n",iter,angle);
            printf("This implies vectors are linearly dependant.\n Restarting the algorithm from current iteration\n");
            fprintf(out_cg,"\nAlgorithm restart at iteration no. %d\n",iter);
            restart_pt=iter;
            iter=0;
            continue;
        }
        else if(iter>0 && angle>5)
        {
            printf("\nAt iteration no. %d , angle between unit direction vectors is %lf degrees\n",iter,angle);
            printf("This implies vectors are linearly independant. Proceed to step 5\n");
        }
        //Step 3(for iteration 0) / Step 5(for other iterations)
        printf("\nPerforming unidirectional search along s(%d) from x(%d)",iter,iter);
        //Calculate bounds for alpha
        double a1,a2,a_lowlim[n],a_uplim[n],a_min,a_max,limits[2];
        for(i=0;i<n;i++)
        {
            if(s[i]==0) // to avoid division by zero situation
            {
                a1=a_lowlim[i-1];
                a2=a_uplim[i-1];
            }
            else
            {
                a1=(lowlim[i]-x[i])/s[i];
                a2=(uplim[i]-x[i])/s[i];
            }

            if(a1<0) //if a1 is negative and a2 is positive then a1 is lower limit and a2 becomes upper limit
            {
                a_lowlim[i]=a1;
                a_uplim[i]=a2;
            }
            else   //if a2 is negative and a1 is positive then a2 is lower limit and a1 becomes upper limit
            {
                a_lowlim[i]=a2;
                a_uplim[i]=a1;
            }

            if(i==0)  // for first iteration the obtained alpha values are stored as a_min and a_max
            {
                a_min=a_lowlim[i];
                a_max=a_uplim[i];
            }

            if(pow(a_lowlim[i],2)<pow(a_min,2)) // if magnitude of current alpha_lower is less than alpha_min update
            {
                a_min=a_lowlim[i];
            }

            if(pow(a_uplim[i],2)<pow(a_max,2)) // if magnitude of current alpha_upper is less than alpha_max update
            {
                a_max=a_uplim[i];
            }
            //printf("\n lowlim %lf \t uplim %lf \t amin %lf \t amax %lf \n",a_lowlim[i],a_uplim[i],a_min,a_max);
        }
        limits[0]=a_min;
        limits[1]=a_max;
        feval2=bounding_phase(x,s,limits,n,R);  // Carrying out unidirectional search by calling bounding phase method
        alpha=newton_raphson(limits[0],limits[1],x,s,feval2,n,R); // Calling newton raphson method
        fevalcg=fevalcg+feval2;
        normx_diffsq=0;
        normx_sq=0;
        printf("\n---------------------------------------------------------------------------------\n");
        printf("\nAfter unidirectional search we get value of x(%d) as : \n",iter+1);
        for(i=0;i<n;i++)
        {
            x_old[i]=x[i];
            normx_sq=normx_sq+(x_old[i]*x_old[i]); //Computing norm of x(k) for checking termination condition
            x[i]=x[i]+(alpha*s[i]);         // Calculating new value of x(k)
            printf("%lf \t",x[i]);
            normx_diffsq=normx_diffsq + pow((x[i]-x_old[i]),2.0);// Computing norm of x(k+1)-x(k)
        }
        //Step 6
        f0=pen_fun(x,n,R,g);
        iter=iter+1;
        fprintf(out_cg,"%d\t %lf\t ", iter,f0);
            for(i=0;i<n;i++)
            {
                fprintf(out_cg,"%lf \t",x[i]);
                fprintf(out_cg,"%lf \t",grad_fx[i]);
            }
            printf("\n Constraint violations : ");
            for(i=0;i<m;i++)
            {
                fprintf(out_cg,"%lf \t",g[i]);
                printf("g%d %lf \t",i+1,g[i]);
            }
            fprintf(out_cg,"%lf\t %lf \t %lf\n",angle,sqrt(normgrad_sqnew),sqrt(normx_diffsq/normx_sq));
            printf("\n\n||x(%d)-x(%d)||/||(x(%d)||\t||grad_f(x(%d))||",iter,iter-1,iter-1,iter);
            printf("\n%lf \t %lf\n",sqrt(normx_diffsq/normx_sq),sqrt(normgrad_sqnew));
            printf("\n Iteration %d is over",iter);
    }while(sqrt(normx_diffsq/normx_sq)>error_2 && sqrt(normgrad_sqnew)>error_3); //Termination Condition
    printf("\n---------------------------------------------------------------------------------\n");
    if(restart_pt>0)
    {
        printf("\nAlgorithm restarted at iteration no. %d\n",restart_pt);
    }
    printf("After %d iterations P(x) = %lf \n", iter,f0);
    printf("The optimum value of x = ( ");
    fprintf(out_cg,"The optimum value of x = ( ");
    for(i=0;i<n;i++)
            {
                fprintf(out_cg,"%lf \t",x[i]);
                printf("%lf \t",x[i]);
            }
    fprintf(out_cg,"\n");
    printf(")\n Constraint violations : ");
    for(i=0;i<m;i++)
            {
                fprintf(out_cg,"%lf \t",g[i]);
                printf("g%d %lf \t",i+1,g[i]);
            }
    printf("\n#Total number of function evaluations: %d",fevalcg);
    fprintf(out_cg,")\n#Total number of function evaluations: %d\n\n",fevalcg);
    fclose(out_cg);
    return f0;
}


