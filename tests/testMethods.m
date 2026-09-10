function tests = testMethods
%TESTMETHODS Numerical tests for the four Krylov methods and their helpers.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repositoryRoot);
testCase.TestData.RepositoryRoot = repositoryRoot;
end

function teardownOnce(testCase)
rmpath(testCase.TestData.RepositoryRoot);
end

function testArnoldiRelation(testCase)
A = gallery("grcar",8);
[V,H,info] = krylov.arnoldi(A,ones(8,1),5);
m = info.Dimension;

verifyFalse(testCase,info.Breakdown);
verifyLessThan(testCase,norm(A*V(:,1:m)-V*H,"fro"),1e-11);
verifyLessThan(testCase,norm(V'*V-eye(size(V,2)),"fro"),1e-11);
end

function testTruncatedArnoldiUsesRecentVectors(testCase)
A = gallery("grcar",10);
k = 2;
[V,AV,info] = krylov.truncatedArnoldi(A,ones(10,1),6,k);

verifyEqual(testCase,AV,A*V,AbsTol=1e-12);
for j = 1:info.Dimension-1
    first = max(1,j-k+1);
    verifyLessThan(testCase,norm(V(:,first:j)'*V(:,j+1)),1e-11);
end
end

function testFunctionHandleOperatorMatchesMatrix(testCase)
A = gallery("grcar",8);
v1 = (1:8)';
[Vm,Hm] = krylov.arnoldi(A,v1,4);
[Vf,Hf] = krylov.arnoldi(@(x) A*x,v1,4);

verifyEqual(testCase,Vf,Vm,AbsTol=1e-13);
verifyEqual(testCase,Hf,Hm,AbsTol=1e-13);
end

function testFullSrdctIsAnIsometry(testCase)
S = krylov.srdct(16,16,2,7);
X = reshape(1:48,16,3);

verifyEqual(testCase,norm(S.Apply(X),"fro"),norm(X,"fro"),RelTol=1e-12);
end

function testSrdctDrawIsDeterministicAndLocal(testCase)
before = rng;
S1 = krylov.srdct(20,9,4,11);
after = rng;
S2 = krylov.srdct(20,9,4,11);

verifyEqual(testCase,after,before);
verifyEqual(testCase,S2.Rows,S1.Rows);
verifyEqual(testCase,S2.Signs,S1.Signs);
verifyEqual(testCase,S2.Apply(ones(20,2)),S1.Apply(ones(20,2)),AbsTol=0);
end

function testSrdctRejectsInvalidSketchSize(testCase)
verifyError(testCase,@() krylov.srdct(10,11,2,0), ...
    "krylov:srdct:InvalidSketchSize");
end

function testPivotedQrRelation(testCase)
B = [eye(5); ones(2,5)];
[Q,R,p] = krylov.pivotedQR(B);

verifyEqual(testCase,B(:,p),Q*R,AbsTol=1e-12);
end

function testGmresFullDimensionMatchesDirectSolve(testCase)
A = gallery("grcar",7) + 2*eye(7);
b = (1:7)';
[x,info] = gmres(A,b,7);

verifyLessThan(testCase,norm(x-A\b)/norm(A\b),1e-11);
verifyLessThan(testCase,info.RelativeResidual,1e-11);
end

function testGmresReturnsExactInitialGuess(testCase)
A = eye(6);
b = ones(6,1);
x0 = b;
[x,info] = gmres(A,b,4,InitialGuess=x0);

verifyEqual(testCase,x,x0,AbsTol=1e-14);
verifyEqual(testCase,info.Dimension,0);
verifyEqual(testCase,info.RelativeResidual,0,AbsTol=1e-14);
end

function testGmresAcceptsFunctionHandle(testCase)
A = gallery("grcar",6) + 2*eye(6);
b = (1:6)';
xMatrix = gmres(A,b,5);
xHandle = gmres(@(v) A*v,b,5);

verifyEqual(testCase,xHandle,xMatrix,AbsTol=1e-12);
end

function testFullSketchSgmresMatchesGmres(testCase)
A = gallery("grcar",10) + 2*eye(10);
b = cos((1:10)');
k = 6;
ell = k;
xClassical = gmres(A,b,k);
[xSketched,info] = sgmres(A,b,k,SketchSize=10,Truncation=ell,Seed=3);

verifyEqual(testCase,xSketched,xClassical,AbsTol=1e-10);
verifyLessThan(testCase,info.SketchedRelativeResidual,1);
end

function testRrFullDimensionMatchesEig(testCase)
A = diag(1:7);
[theta,U,info] = rr(A,ones(7,1),7,NumEigenpairs=3);

verifyEqual(testCase,theta,[7;6;5],AbsTol=1e-11);
verifyLessThan(testCase,max(info.RitzResiduals),1e-11);
verifyEqual(testCase,vecnorm(U),ones(1,3),AbsTol=1e-12);
end

function testFullSketchSrrMatchesRr(testCase)
A = diag(1:8);
v1 = (1:8)';
k = 6;
ell = k;
[thetaClassical] = rr(A,v1,k,NumEigenpairs=3);
[thetaSketched,~,info] = srr(A,v1,k,SketchSize=8,Truncation=ell, ...
    NumEigenpairs=3,Seed=5);

verifyEqual(testCase,thetaSketched,thetaClassical,AbsTol=1e-10);
verifyLessThan(testCase,max(info.SketchedRitzResiduals),1);
end

function testMethodsReportNestedTiming(testCase)
A = gallery("grcar",8)+2*eye(8);
b = (1:8)';

[~,gmresInfo] = gmres(A,b,5);
[~,sgmresInfo] = sgmres(A,b,5,Seed=2);
[~,~,rrInfo] = rr(A,b,5,NumEigenpairs=1);
[~,~,srrInfo] = srr(A,b,5,NumEigenpairs=1,Seed=2);
[~,exactInfo] = gmres(eye(8),ones(8,1),5, ...
    InitialGuess=ones(8,1));

infos = {gmresInfo,sgmresInfo,rrInfo,srrInfo,exactInfo};
for index = 1:numel(infos)
    info = infos{index};
    verifyTrue(testCase,isfinite(info.CoreTime));
    verifyTrue(testCase,isfinite(info.WallTime));
    verifyGreaterThanOrEqual(testCase,info.CoreTime,0);
    verifyGreaterThanOrEqual(testCase,info.WallTime,info.CoreTime);
end
end
