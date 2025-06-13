#' GhostBASIL solver
#'
#' This function fits a lasso regression model between disease of interest and gene expression levels (observed and knockoff copies). It solves the quadratic optimization problem
#' \deqn{\min_{\beta}\frac{1}{2}\cdot\beta^TG\beta-\beta^Tb+\lambda\|\beta\|_1,}
#' over a set of tuning parameters \eqn{\lambda}'s, where
#' \eqn{G=(1_{M+1}1_{M+1}^T)\otimes\Sigma_{S}+I_{M+1}\otimes S}.
#'
#' @param b         Coefficient of linear term.
#' @param Sigma_S   Correlation matrix \eqn{Sigma} minus the diagonal matrix \eqn{S}
#' @param S         Diagnoal elements of the S-matrix.
#' @param M         The number of knockoff copies (under the multiple knockoff framework).
#' @param lambda    User-sepcific tuning parameter list.
#'
#' @return A list containing all considered tuning parameters \eqn{\lambda}'s and the corresponing solutions.
#' @importFrom Matrix Matrix
#' @importFrom Matrix colSums
#' @importFrom Matrix t
#' @importFrom Rfast spdinv
#' @importFrom methods as
#' @importFrom methods cbind2
#' @importFrom methods rbind2
#' @export
Ghost_Basil<-function(b,Sigma_S,S,M,lambda=NULL)
{
  dim_all<-c(length(b),dim(Sigma_S)*(M+1),length(S)*(M+1))
  if(max(dim_all)!=min(dim_all))
  {
    print("Error: Dimension inconsistency!")
  }
  else
  {
    p<-length(S)

    gradient_f<--b
    lambda_max<-max(abs(gradient_f))
    lambda_min<-lambda_max/1000
    lambda_seq<-seq(lambda_max,lambda_min,(lambda_min-lambda_max)/1000)
    if(!is.null(lambda))
    {
      lambda_seq<-sort(unique(c(lambda_seq[lambda_seq>=min(lambda)],lambda)),decreasing = T)
    }

    solution_path<-Matrix(0,nrow = (M+1)*p,ncol = sum(lambda_seq>=lambda_max))
    solution_current<-solution_path[,ncol(solution_path)]

    S_zero<-1:((M+1)*p)
    S_nonzero<-vector()
    S_new<-which.max(abs(gradient_f))
    S_zero<-setdiff(S_zero,S_new)
    S_nonzero<-c(S_nonzero,S_new)
    S_deactive<-NULL
    Span_deactive<-NULL

    delta_current<-rep(0,length(b))
    delta_current[S_new]<--sign(gradient_f[S_new])
    count<-ncol(solution_path)

    remainders_nonzero<-S_new%%p
    L_nonzero_inv<-matrix(1/sqrt(Sigma_S[remainders_nonzero,remainders_nonzero]+S[remainders_nonzero]),nrow = 1,ncol=1)
    L_nonzero_inv<-as(L_nonzero_inv, "dgCMatrix")
    L_nonzero_inv_b_nonzero<-as.vector(L_nonzero_inv*b[S_new])
    L_nonzero_inv_sign<-as.vector(L_nonzero_inv*delta_current[S_new])

    while ((count<length(lambda_seq)))
      # while ((count<838))
    {
      count<-1+count

      solution_nonzero_prop<-colSums(L_nonzero_inv*(L_nonzero_inv_b_nonzero-lambda_seq[count]*L_nonzero_inv_sign))
      solution_nonzero<-solution_nonzero_prop
      if(length(S_deactive)!=0)
      {
        Span_deactive<-t(L_nonzero_inv)%*%L_nonzero_inv[,S_deactive]
        if(length(S_deactive)==1)
        {adjust_factor<-Span_deactive%*%(solve(Span_deactive[S_deactive,])%*%solution_nonzero_prop[S_deactive])}
        else
        {adjust_factor<-Span_deactive%*%(spdinv(as.matrix(Span_deactive[S_deactive,]))%*%solution_nonzero_prop[S_deactive])}

        solution_nonzero<-solution_nonzero_prop-adjust_factor[,1]
      }

      # S_deactive<-NULL
      S_deactive_new<-setdiff(which(sign(solution_nonzero)!=delta_current[S_nonzero]),S_deactive)

      while(length(S_deactive_new)!=0)
      {

        S_deactive<-c(S_deactive,S_deactive_new)

        Span_deactive<-t(L_nonzero_inv)%*%L_nonzero_inv[,S_deactive]
        if(length(S_deactive)==1)
        {adjust_factor<-Span_deactive%*%(solve(Span_deactive[S_deactive,])%*%solution_nonzero_prop[S_deactive])}
        else
        {adjust_factor<-Span_deactive%*%(spdinv(as.matrix(Span_deactive[S_deactive,]))%*%solution_nonzero_prop[S_deactive])}

        solution_nonzero<-solution_nonzero_prop-adjust_factor[,1]

        L_nonzero_inv_sign<-L_nonzero_inv_sign-colSums(t(L_nonzero_inv[,S_deactive])*delta_current[S_nonzero[S_deactive]])
        delta_current[S_nonzero[S_deactive]]<-0
        solution_nonzero[S_deactive]<-0

        S_deactive_new<-which(sign(solution_nonzero)!=delta_current[S_nonzero])
      }

      solution_new<-solution_current
      solution_new[S_nonzero]<-solution_nonzero

      #Checkpoint
      gradient_f_new<-rep(colSums(t(Sigma_S[,remainders_nonzero])*solution_new[S_nonzero]),M+1)+rep(S,M+1)*solution_new-b

      if(max(abs(gradient_f_new[c(S_zero)]))<=lambda_seq[count])
      {
        solution_path<-cbind(solution_path,solution_new)
        gradient_f<-gradient_f_new
        solution_current<-solution_new
        # gradient_f_all<-cbind(gradient_f_all,solution_new)
      }
      else
      {
        # break
        count<-count-1
        {
          S_new<-S_zero[which.max(abs(gradient_f_new[S_zero]))]
          remainders_new<-S_new%%p

          bb<-Sigma_S[remainders_nonzero,remainders_new]
          cc<-Sigma_S[remainders_new,remainders_new]+S[remainders_new]

          L_b<-(L_nonzero_inv%*%bb)[,1]
          cc_L_b<-sqrt(cc-sum(L_b^2))
          newcolumn_L_nonzero_inv<-c(-colSums(L_nonzero_inv*L_b)/cc_L_b,1/cc_L_b)


          delta_current[S_new]<--sign(gradient_f[S_new])
          S_zero<-setdiff(S_zero,S_new)
          S_nonzero<-c(S_nonzero,S_new)
          remainders_nonzero<-c(remainders_nonzero,remainders_new)

          L_nonzero_inv<-rbind2(cbind2(L_nonzero_inv,0),newcolumn_L_nonzero_inv)
          L_nonzero_inv_b_nonzero<-c(L_nonzero_inv_b_nonzero,sum(newcolumn_L_nonzero_inv*b[S_nonzero]))
          L_nonzero_inv_sign<-c(L_nonzero_inv_sign,sum(newcolumn_L_nonzero_inv*delta_current[S_nonzero]))
        }
      }

      {
        # # print(c(max(abs(gradient_f_new[c(S_zero,S_nonzero[S_deactive])])),lambda_seq[count]))
        # if(max(abs(gradient_f_new[c(S_zero,S_nonzero[S_deactive])]))<=lambda_seq[count])
        # {
        #   solution_path<-cbind(solution_path,solution_new)
        #   gradient_f<-gradient_f_new
        #   solution_current<-solution_new
        #   # gradient_f_all<-cbind(gradient_f_all,solution_new)
        # }else
        # {
        #   # break
        #   count<-count-1
        #   if(length(S_deactive)==0)
        #   {
        #     S_new<-S_zero[which.max(abs(gradient_f_new[S_zero]))]
        #     remainders_new<-S_new%%p
        #
        #     bb<-Sigma_S[remainders_nonzero,remainders_new]
        #     cc<-Sigma_S[remainders_new,remainders_new]+S[remainders_new]
        #
        #     L_b<-(L_nonzero_inv%*%bb)[,1]
        #     cc_L_b<-sqrt(cc-sum(L_b^2))
        #     newcolumn_L_nonzero_inv<-c(-colSums(L_nonzero_inv*L_b)/cc_L_b,1/cc_L_b)
        #
        #
        #     delta_current[S_new]<--sign(gradient_f[S_new])
        #     S_zero<-setdiff(S_zero,S_new)
        #     S_nonzero<-c(S_nonzero,S_new)
        #     remainders_nonzero<-c(remainders_nonzero,remainders_new)
        #
        #     L_nonzero_inv<-rbind2(cbind2(L_nonzero_inv,0),newcolumn_L_nonzero_inv)
        #     L_nonzero_inv_b_nonzero<-c(L_nonzero_inv_b_nonzero,sum(newcolumn_L_nonzero_inv*b[S_nonzero]))
        #     L_nonzero_inv_sign<-c(L_nonzero_inv_sign,sum(newcolumn_L_nonzero_inv*delta_current[S_nonzero]))
        #   }else
        #   {
        #     if(max(abs(gradient_f_new[c(S_nonzero[S_deactive])]))<=lambda_seq[count])
        #     {
        #       count<-count+1
        #       reactive_index<-which.max(abs(gradient_f_new[c(S_nonzero[S_deactive])]))
        #
        #       delta_current[S_nonzero[S_deactive[reactive_index]]]<--sign(gradient_f_new[c(S_nonzero[S_deactive[reactive_index]])])
        #       L_nonzero_inv_sign<-L_nonzero_inv_sign+colSums(t(L_nonzero_inv[,S_deactive[reactive_index]])*delta_current[S_nonzero[S_deactive[reactive_index]]])
        #       S_deactive<-S_deactive[-reactive_index]
        #
        #       solution_path<-cbind(solution_path,solution_new)
        #       gradient_f<-gradient_f_new
        #       solution_current<-solution_new
        #
        #     }else
        #     {
        #       S_new<-S_zero[which.max(abs(gradient_f_new[S_zero]))]
        #       remainders_new<-S_new%%p
        #
        #       bb<-Sigma_S[remainders_nonzero,remainders_new]
        #       cc<-Sigma_S[remainders_new,remainders_new]+S[remainders_new]
        #
        #       L_b<-(L_nonzero_inv%*%bb)[,1]
        #       cc_L_b<-sqrt(cc-sum(L_b^2))
        #       newcolumn_L_nonzero_inv<-c(-colSums(L_nonzero_inv*L_b)/cc_L_b,1/cc_L_b)
        #
        #
        #       delta_current[S_new]<--sign(gradient_f[S_new])
        #       S_zero<-setdiff(S_zero,S_new)
        #       S_nonzero<-c(S_nonzero,S_new)
        #       remainders_nonzero<-c(remainders_nonzero,remainders_new)
        #
        #       L_nonzero_inv<-rbind2(cbind2(L_nonzero_inv,0),newcolumn_L_nonzero_inv)
        #       L_nonzero_inv_b_nonzero<-c(L_nonzero_inv_b_nonzero,sum(newcolumn_L_nonzero_inv*b[S_nonzero]))
        #       L_nonzero_inv_sign<-c(L_nonzero_inv_sign,sum(newcolumn_L_nonzero_inv*delta_current[S_nonzero]))
        #     }
        #   }
        # }
      }

      # print(c(count,length(S_nonzero)))
    }
    if(!is.null(lambda))
    {
      # gamma_seq<-sort(c(gamma_seq,gamma),decreasing = T)
      solution_path<-solution_path[,match(lambda,lambda_seq)]
      lambda_seq<-lambda_seq[match(lambda,lambda_seq)]
    }
    return(list(lambda=lambda_seq,Solution=solution_path))
  }
}
