// ClothMesh.cpp

#include "ClothMesh.h"
#include "ProceduralMeshComponent.h"
#include "Engine/World.h"
#include "GameFramework/PlayerController.h"
#include "Components/InputComponent.h"
#include "DrawDebugHelpers.h"
#include "Kismet/KismetMathLibrary.h"
#include "Engine/Engine.h"
#include "Engine/CollisionProfile.h"
#include "Math/UnrealMathUtility.h"

// Sets default values
AClothMesh::AClothMesh()
{
    PrimaryActorTick.bCanEverTick = true;

    // Initialize the Procedural Mesh Component
    ClothMesh = CreateDefaultSubobject<UProceduralMeshComponent>(TEXT("ClothMesh"));
    ClothMesh->bUseAsyncCooking = false;
    ClothMesh->SetCollisionEnabled(ECollisionEnabled::QueryAndPhysics);
    ClothMesh->SetCollisionObjectType(ECC_WorldDynamic);
    ClothMesh->SetCollisionResponseToAllChannels(ECR_Block);
    ClothMesh->SetSimulatePhysics(false);
    ClothMesh->bCastDynamicShadow = true;
    ClothMesh->SetGenerateOverlapEvents(true);

    // Set root component
    SetRootComponent(ClothMesh);

    // Initialize variables
    ClothWidth = 30;
    ClothHeight = 30;
    ClothSpacing = 20.0f;
    NumberOfHooks = 6; // Ensure this is even
    if (NumberOfHooks % 2 != 0)
    {
        NumberOfHooks += 1; // Make it even if it's odd
    }
    GravityStrength = -980.0f;
    WindStrength = 900;
    WindDirection = FVector(1, 0, 0); // Wind blowing along the X-axis
    Damping = 0.99f;
    SelfCollisionRadius = 15.0f;
    SelfCollisionStiffness = 0.7f;
    TearThreshold = 1.5f;
    bTearingEnabled = false;

    MinHookSpacing = 50.0f;

    // Initialize curtain variables
    bMovingLeft = false;
    bMovingRight = false;
    bIsOpening = false;
    bIsClosing = false;
    PairMovementSpeed = 200.0f;

    CurrentClosingPairIndex = 0;
    CurrentOpeningPairIndex = 0;
    TotalPairs = NumberOfHooks / 2;

    CurrentOpeningHookIndex = 0;
    CurrentClosingHookIndex = 0;
    bIsMovingHook = false;
}

// Called when the game starts or when spawned
void AClothMesh::BeginPlay()
{
    Super::BeginPlay();

    InitializeParticlesAndSprings();

    // Setup input bindings
    EnableInput(GetWorld()->GetFirstPlayerController());
    if (InputComponent)
    {
        InputComponent->BindAction("ResetCloth", IE_Pressed, this, &AClothMesh::ResetCloth);
        InputComponent->BindAction("GrabCloth", IE_Pressed, this, &AClothMesh::OnMouseClick);
        // Bind curtain controls if needed
    }

    // Initialize curtain targets
    InitializeCurtainTargets();
}

// Called every frame
void AClothMesh::Tick(float DeltaTime)
{
    Super::Tick(DeltaTime);

    // Lock top rows to prevent Z movement
    LockTopRowZAxis();

    // Simulate the cloth and update the mesh based on particle positions
    SimulateCloth(DeltaTime);
    UpdateMesh();

    // Update curtain movement if in progress
    UpdateCurtainMovement(DeltaTime);

    // Debug: Draw spheres at each hook position to represent pinned points
    for (int32 PinnedIndex : PinnedParticles)
    {
        DrawDebugSphere(
            GetWorld(),
            Particles[PinnedIndex].Position, // Position of the hook
            10.0f,                           // Radius of the sphere
            12,                              // Number of segments
            FColor::Green,                   // Color of the sphere
            false,                           // Persistent lines (false = only for this frame)
            -1.0f,                           // Lifetime (negative = only for this frame)
            0,                               // Depth priority
            2.0f                             // Line thickness
        );
    }
}

void AClothMesh::InitializeParticlesAndSprings()
{
    // Clear existing arrays
    Particles.Empty();
    Springs.Empty();
    PinnedParticles.Empty();
    InitialHookPositions.Empty();
    InitialHookPositionsCurtain.Empty();

    // Initialize particles above ground
    float InitialZ = 1000.0f; // 10 meters above ground (assuming 1 unit = 1 cm)

    // Ensure a minimum of 2 hooks for pairing
    int32 HooksToCreate = FMath::Max(NumberOfHooks, 2);
    int32 MidIndex = ClothWidth / 2;

    // Determine indices for hooks: left end, right end, and additional evenly spaced hooks
    TArray<int32> HookIndices;

    // Always add leftmost and rightmost hooks
    HookIndices.Add(0); // Leftmost
    HookIndices.Add(ClothWidth - 1); // Rightmost

    int32 HooksToAdd = HooksToCreate - 2;

    // Evenly distribute the remaining hooks between left and right
    for (int32 i = 1; i <= HooksToAdd; i++)
    {
        float fraction = static_cast<float>(i) / (HooksToAdd + 1);
        int32 X = FMath::RoundToInt(fraction * (ClothWidth - 1));
        HookIndices.Add(X);
    }

    // Ensure HookIndices has exactly NumberOfHooks elements
    HookIndices.SetNum(NumberOfHooks);

    // Remove duplicates and sort HookIndices
    HookIndices.Sort();
    for (int32 i = 0; i < HookIndices.Num(); ++i)
    {
        int32 CurrentValue = HookIndices[i];

        // Loop backwards from the end of the array to the element after the current one
        for (int32 j = HookIndices.Num() - 1; j > i; --j)
        {
            if (HookIndices[j] == CurrentValue)
            {
                // Remove duplicate found at index j
                HookIndices.RemoveAt(j);
            }
        }
    }

    // Log HookIndices for debugging
    UE_LOG(LogTemp, Warning, TEXT("HookIndices:"));
    for (auto Hook : HookIndices)
    {
        UE_LOG(LogTemp, Warning, TEXT("%d"), Hook);
    }

    // Ensure at least 2 hooks
    if (HookIndices.Num() < 2)
    {
        UE_LOG(LogTemp, Warning, TEXT("Not enough hooks. Adding hooks to meet the minimum requirement."));
        // Add hooks to meet the minimum of 2
        if (!HookIndices.Contains(0))
            HookIndices.Add(0);
        if (!HookIndices.Contains(ClothWidth - 1))
            HookIndices.Add(ClothWidth - 1);
    }

    // Create particles
    for (int32 Y = 0; Y < ClothHeight; Y++)
    {
        for (int32 X = 0; X < ClothWidth; X++)
        {
            FParticle Particle;
            Particle.Position = FVector(X * ClothSpacing, Y * ClothSpacing, InitialZ);
            Particle.PreviousPosition = Particle.Position;
            Particle.Acceleration = FVector::ZeroVector;

            // Pin particles based on HookIndices
            if (Y == 0 && HookIndices.Contains(X))
            {
                Particle.bIsPinned = true;
                PinnedParticles.Add(Y * ClothWidth + X); // Index of the pinned particle

                // Store the initial position for each pinned particle (hook)
                InitialHookPositions.Add(Particle.Position);
                InitialHookPositionsCurtain.Add(Particle.Position); // For curtain reset
            }
            else
            {
                Particle.bIsPinned = false;
            }

            Particles.Add(Particle);
        }
    }

    // Sort PinnedParticles based on their X position (left to right)
    PinnedParticles.Sort([this](const int32 A, const int32 B) -> bool {
        return Particles[A].Position.X < Particles[B].Position.X;
        });

    // Store the sorted initial positions
    InitialHookPositions.Empty();
    InitialHookPositionsCurtain.Empty();
    for (int32 PinnedIndex : PinnedParticles)
    {
        InitialHookPositions.Add(Particles[PinnedIndex].Position);
        InitialHookPositionsCurtain.Add(Particles[PinnedIndex].Position); // For curtain reset
    }

    // Initialize curtain target positions
    InitializeCurtainTargets();

    // Create springs (structural, shear, and bend)
    for (int32 Y = 0; Y < ClothHeight; Y++)
    {
        for (int32 X = 0; X < ClothWidth; X++)
        {
            int32 Index = X + Y * ClothWidth;

            /// Structural springs
            if (X < ClothWidth - 1)
                Springs.Add(FSpring(Index, Index + 1, ClothSpacing, 1.5f)); 
            if (Y < ClothHeight - 1)
                Springs.Add(FSpring(Index, Index + ClothWidth, ClothSpacing, 1.5f)); 

            // Shear springs
            if (X < ClothWidth - 1 && Y < ClothHeight - 1)
            {
                Springs.Add(FSpring(Index, Index + ClothWidth + 1, FMath::Sqrt(2) * ClothSpacing, 1.0f)); 
                Springs.Add(FSpring(Index + 1, Index + ClothWidth, FMath::Sqrt(2) * ClothSpacing, 1.0f));
            }

            // Bend springs
            if (X < ClothWidth - 2)
                Springs.Add(FSpring(Index, Index + 2, 2 * ClothSpacing, 0.8f));
            if (Y < ClothHeight - 2)
                Springs.Add(FSpring(Index, Index + 2 * ClothWidth, 2 * ClothSpacing, 0.8f)); 


            // Additional springs for tightly woven structure if enabled
            if (ClothWidth > 3 && ClothHeight > 3)
            {
                if (X < ClothWidth - 2 && Y < ClothHeight - 2)
                {
                    Springs.Add(FSpring(Index, Index + ClothWidth + 2, FMath::Sqrt(5) * ClothSpacing, 1.2f)); // Adjusted
                    Springs.Add(FSpring(Index + 2, Index + ClothWidth, FMath::Sqrt(5) * ClothSpacing, 1.2f)); // Adjusted
                }
            }
        }
    }

    // Initialize movement indices and target positions
    CurrentOpeningPairIndex = 0;
    CurrentClosingPairIndex = 0;
    TargetHookPositionsLeft.Empty();
    TargetHookPositionsRight.Empty();

    // Generate the initial mesh for the cloth
    GenerateClothMesh();
}

void AClothMesh::GenerateClothMesh()
{
    VertexPositions.Empty();
    TArray<int32> Triangles;
    TArray<FVector> Normals;
    TArray<FVector2D> UVs;
    TArray<FLinearColor> VertexColors;
    TArray<FProcMeshTangent> Tangents;

    for (const FParticle& Particle : Particles)
    {
        VertexPositions.Add(Particle.Position);
        UVs.Add(FVector2D((Particle.Position.X) / (ClothWidth * ClothSpacing), (Particle.Position.Y) / (ClothHeight * ClothSpacing)));
        Normals.Add(FVector(0, 0, 1));
        VertexColors.Add(FLinearColor::White);
        Tangents.Add(FProcMeshTangent(1, 0, 0));
    }

    // Generate the triangles
    for (int32 Y = 0; Y < ClothHeight - 1; Y++)
    {
        for (int32 X = 0; X < ClothWidth - 1; X++)
        {
            int32 Index = X + Y * ClothWidth;

            // First triangle of the quad
            Triangles.Add(Index);
            Triangles.Add(Index + ClothWidth);
            Triangles.Add(Index + 1);

            // Second triangle of the quad
            Triangles.Add(Index + 1);
            Triangles.Add(Index + ClothWidth);
            Triangles.Add(Index + ClothWidth + 1);
        }
    }

    // Create the mesh section
    ClothMesh->CreateMeshSection_LinearColor(0, VertexPositions, Triangles, Normals, UVs, VertexColors, Tangents, true);
    ClothMesh->SetCollisionEnabled(ECollisionEnabled::QueryAndPhysics);
    ClothMesh->SetSimulatePhysics(false);
}

void AClothMesh::SimulateCloth(float DeltaTime)
{
    TimeStep = DeltaTime;

    // Apply external forces like gravity and wind
    ApplyForces(DeltaTime);

    // Verlet Integration with velocity tracking
    for (int32 i = 0; i < Particles.Num(); i++)
    {
        if (Particles[i].bIsPinned)
            continue; // Skip updating position for pinned particles

        FVector Temp = Particles[i].Position;
        // Calculate velocity based on current and previous positions
        Particles[i].Velocity = (Particles[i].Position - Particles[i].PreviousPosition) / DeltaTime * Damping;
        // Update position with velocity and acceleration
        Particles[i].Position += Particles[i].Velocity * DeltaTime + Particles[i].Acceleration * DeltaTime * DeltaTime;

        // Limit maximum velocity to prevent instability
        float MaxVelocity = 500.0f;
        if (Particles[i].Velocity.Size() > MaxVelocity)
        {
            Particles[i].Velocity = Particles[i].Velocity.GetSafeNormal() * MaxVelocity;
            Particles[i].Position = Particles[i].PreviousPosition + Particles[i].Velocity * DeltaTime;
        }

        Particles[i].PreviousPosition = Temp;
    }

    // Enforce Fixed Z Positions for Hooks More Strictly
    for (int32 i = 0; i < PinnedParticles.Num(); i++)
    {
        int32 Index = PinnedParticles[i];
        // Completely lock the position except for X and Y
        Particles[Index].Position = FVector(Particles[Index].Position.X, Particles[Index].Position.Y, InitialHookPositions[i].Z);
        Particles[Index].PreviousPosition = FVector(Particles[Index].PreviousPosition.X, Particles[Index].PreviousPosition.Y, InitialHookPositions[i].Z);
    }

    // Resolve constraints multiple times for stability
    for (int32 Iter = 0; Iter < 10; Iter++) // Increased iterations for better stability
    {
        // Resolve springs
        for (int32 s = 0; s < Springs.Num(); s++)
        {
            FSpring& Spring = Springs[s];
            FParticle& pA = Particles[Spring.ParticleA];
            FParticle& pB = Particles[Spring.ParticleB];

            FVector Delta = pB.Position - pA.Position;
            float CurrentLength = Delta.Size();
            if (CurrentLength == 0.0f)
                continue; // Prevent division by zero

            float Difference = (CurrentLength - Spring.RestLength) / CurrentLength;

            if (pA.bIsPinned && pB.bIsPinned)
                continue;

            if (!pA.bIsPinned && !pB.bIsPinned)
            {
                pA.Position += Delta * 0.5f * Spring.Stiffness * Difference;
                pB.Position -= Delta * 0.5f * Spring.Stiffness * Difference;
            }
            else if (pA.bIsPinned)
            {
                pB.Position -= Delta * Spring.Stiffness * Difference;
            }
            else if (pB.bIsPinned)
            {
                pA.Position += Delta * Spring.Stiffness * Difference;
            }
        }

        // Handle collisions with ground and other objects
        HandleCollisions();

        // Handle self-collisions within the cloth
        HandleSelfCollisions();

        // Handle tearing if enabled
        if (bTearingEnabled)
        {
            HandleTearing();
        }
    }
}

void AClothMesh::HandleCollisions()
{
    float GroundLevel = 0.0f; // Ground is at Z = 0
    const float FrictionCoefficient = 0.8f; // Adjust as needed
    const float Restitution = 0.3f; // Bounciness

    for (FParticle& Particle : Particles)
    {
        if (Particle.bIsPinned)
            continue; // Skip collision handling for pinned particles

        // Ground collision
        if (Particle.Position.Z < GroundLevel)
        {
            // Calculate penetration depth
            float Penetration = GroundLevel - Particle.Position.Z;

            // Correct the position based on penetration
            Particle.Position.Z += Penetration;

            // Reflect the Z component of velocity based on restitution
            FVector Velocity = (Particle.Position - Particle.PreviousPosition) / TimeStep;
            Velocity.Z = -Velocity.Z * Restitution;
            Particle.PreviousPosition = Particle.Position - Velocity * TimeStep;

            // Apply friction to X and Y velocities
            FVector HorizontalVelocity = FVector(Velocity.X, Velocity.Y, 0.0f);
            HorizontalVelocity *= FrictionCoefficient;
            Particle.PreviousPosition.X = Particle.Position.X - HorizontalVelocity.X * TimeStep;
            Particle.PreviousPosition.Y = Particle.Position.Y - HorizontalVelocity.Y * TimeStep;
        }

        // Collision with other objects (if any)
        for (AActor* CollisionActor : CollisionObjects)
        {
            if (!CollisionActor)
                continue;

            UPrimitiveComponent* PrimComp = CollisionActor->FindComponentByClass<UPrimitiveComponent>();
            if (!PrimComp)
                continue;

            FVector ClosestPoint;
            float Distance = PrimComp->GetClosestPointOnCollision(Particle.Position, ClosestPoint, NAME_None);
            FVector Delta = Particle.Position - ClosestPoint;

            if (Distance < SelfCollisionRadius && Distance > 0.0f) // Using SelfCollisionRadius as a generic collision threshold
            {
                FVector CollisionNormal = Delta.GetSafeNormal();
                FVector Correction = CollisionNormal * (SelfCollisionRadius - Distance) * SelfCollisionStiffness;

                // Position correction
                if (!Particle.bIsPinned)
                {
                    Particle.Position += Correction;
                    Particle.PreviousPosition += Correction;
                }
            }
        }
    }
}

void AClothMesh::ApplyForces(float DeltaTime)
{
    // Apply gravity
    FVector Gravity = FVector(0, 0, GravityStrength); 

    // Apply constant wind
    FVector Wind = WindDirection.GetSafeNormal() * WindStrength;

    for (FParticle& Particle : Particles)
    {
        if (Particle.bIsPinned) 
            continue; // Skip applying forces to pinned particles

        // Reset acceleration
        Particle.Acceleration = FVector::ZeroVector;

        // Apply gravity
        Particle.Acceleration += Gravity;

        // Apply wind
        Particle.Acceleration += Wind;
    }
}

int32 GetCellHash(const FVector& Position, float CellSize)
{
    int32 x = FMath::FloorToInt(Position.X / CellSize);
    int32 y = FMath::FloorToInt(Position.Y / CellSize);
    int32 z = FMath::FloorToInt(Position.Z / CellSize);
    return x + y * 73856093 + z * 19349663; // Large prime numbers for hashing
}

void AClothMesh::HandleSelfCollisions()
{
    float CellSize = SelfCollisionRadius * 2.0f;
    TMap<FIntVector, TArray<int32>> SpatialHashVector;

    // Populate the spatial hash with cell coordinates
    for (int32 i = 0; i < Particles.Num(); i++)
    {
        FVector Pos = Particles[i].Position;
        FIntVector Cell = FIntVector(
            FMath::FloorToInt(Pos.X / CellSize),
            FMath::FloorToInt(Pos.Y / CellSize),
            FMath::FloorToInt(Pos.Z / CellSize)
        );
        SpatialHashVector.FindOrAdd(Cell).Add(i);
    }

    // Define neighbor offsets for a 3D grid (including the cell itself)
    TArray<FIntVector> NeighborOffsets;
    for (int32 x = -1; x <= 1; x++)
    {
        for (int32 y = -1; y <= 1; y++)
        {
            for (int32 z = -1; z <= 1; z++)
            {
                NeighborOffsets.Add(FIntVector(x, y, z));
            }
        }
    }

    // Iterate through each cell and its neighbors
    for (const auto& Elem : SpatialHashVector)
    {
        FIntVector Cell = Elem.Key;
        const TArray<int32>& CellParticles = Elem.Value;

        for (const FIntVector& Offset : NeighborOffsets)
        {
            FIntVector NeighborCell = Cell + Offset;

            // Check if neighbor cell exists
            if (!SpatialHashVector.Contains(NeighborCell))
                continue;

            const TArray<int32>& NeighborParticles = SpatialHashVector[NeighborCell];

            for (int32 i = 0; i < CellParticles.Num(); i++)
            {
                int32 IndexA = CellParticles[i];
                FParticle& pA = Particles[IndexA];

                for (int32 j = 0; j < NeighborParticles.Num(); j++)
                {
                    int32 IndexB = NeighborParticles[j];

                    // Prevent duplicate checks and self-collision
                    if (IndexA >= IndexB)
                        continue;

                    FParticle& pB = Particles[IndexB];

                    FVector Delta = pB.Position - pA.Position;
                    float Distance = Delta.Size();

                    if (Distance < SelfCollisionRadius && Distance > 0.0f)
                    {
                        FVector CollisionNormal = Delta.GetSafeNormal();
                        FVector RelativeVelocity = pB.Velocity - pA.Velocity;

                        // Compute relative velocity components
                        float Vn = FVector::DotProduct(RelativeVelocity, CollisionNormal);
                        if (Vn > 0)
                            continue; // Prevent particles moving away from each other

                        FVector VnVec = Vn * CollisionNormal;
                        FVector VtVec = RelativeVelocity - VnVec;

                        // Apply friction to tangential component
                        FVector VtDir = VtVec.GetSafeNormal();
                        FVector FrictionForce = -VtDir * DynamicFrictionCoefficient;

                        // Limit the friction force to not reverse the tangential velocity
                        if (FVector::DotProduct(FrictionForce, VtVec) > 0.0f)
                        {
                            FrictionForce = -VtVec;
                        }

                        // Apply friction to velocities
                        pA.Velocity += FrictionForce * TimeStep;
                        pB.Velocity -= FrictionForce * TimeStep;

                        // Position correction to prevent overlap
                        FVector Correction = CollisionNormal * (SelfCollisionRadius - Distance) * 0.5f * SelfCollisionStiffness;

                        if (!pA.bIsPinned && !pB.bIsPinned)
                        {
                            pA.Position -= Correction;
                            pB.Position += Correction;
                        }
                        else if (pA.bIsPinned && !pB.bIsPinned)
                        {
                            pB.Position += Correction * 2.0f; // Apply more correction to unpinned
                        }
                        else if (!pA.bIsPinned && pB.bIsPinned)
                        {
                            pA.Position -= Correction * 2.0f; // Apply more correction to unpinned
                        }
                        // If both are pinned, no correction is applied
                    }
                }
            }
        }
    }
}

void AClothMesh::UpdateMesh()
{
    // Update vertex positions based on particle positions
    for (int32 i = 0; i < Particles.Num(); i++)
    {
        VertexPositions[i] = Particles[i].Position;
    }

    // Recalculate normals for smooth shading
    TArray<FVector> Normals;
    Normals.Init(FVector::ZeroVector, VertexPositions.Num());

    // Iterate through each triangle to calculate normals
    for (int32 Y = 0; Y < ClothHeight - 1; Y++)
    {
        for (int32 X = 0; X < ClothWidth - 1; X++)
        {
            int32 Index = X + Y * ClothWidth;

            // First triangle of the quad
            FVector v0 = VertexPositions[Index];
            FVector v1 = VertexPositions[Index + ClothWidth];
            FVector v2 = VertexPositions[Index + 1];

            FVector Edge1 = v1 - v0;
            FVector Edge2 = v2 - v0;
            FVector Normal1 = FVector::CrossProduct(Edge1, Edge2).GetSafeNormal();

            Normals[Index] += Normal1;
            Normals[Index + ClothWidth] += Normal1;
            Normals[Index + 1] += Normal1;

            // Second triangle of the quad
            FVector v3 = VertexPositions[Index + ClothWidth + 1];

            FVector Edge3 = v2 - v1;
            FVector Edge4 = v3 - v1;
            FVector Normal2 = FVector::CrossProduct(Edge3, Edge4).GetSafeNormal();

            Normals[Index + 1] += Normal2;
            Normals[Index + ClothWidth] += Normal2;
            Normals[Index + ClothWidth + 1] += Normal2;
        }
    }

    // Normalize the normals
    for (FVector& Normal : Normals)
    {
        Normal = Normal.GetSafeNormal();
    }

    // Update the mesh section with new positions and normals
    ClothMesh->UpdateMeshSection(0, VertexPositions, Normals, {}, {}, {});
}

void AClothMesh::ResetCloth()
{
    InitializeParticlesAndSprings();
    CurrentOpeningPairIndex = 0;
    CurrentClosingPairIndex = 0;
    TargetHookPositionsLeft.Empty();
    TargetHookPositionsRight.Empty();
    bMovingLeft = false;
    bMovingRight = false;
    bIsOpening = false;
    bIsClosing = false;

    // Reinitialize curtain targets
    InitializeCurtainTargets();
}

void AClothMesh::ReleaseCloth()
{
    for (int32 PinnedIndex : PinnedParticles)
    {
        Particles[PinnedIndex].bIsPinned = false;
    }
    PinnedParticles.Empty(); // Clear the array as hooks are no longer pinned
}

void AClothMesh::SetWind(FVector NewWindDirection, float NewWindStrength)
{
    WindDirection = NewWindDirection;
    WindStrength = NewWindStrength;
}

void AClothMesh::ToggleTightlyWoven(bool bEnable)
{
    if (bEnable)
    {
        // Add additional springs for tighter weaving
        for (int32 Y = 0; Y < ClothHeight; Y++)
        {
            for (int32 X = 0; X < ClothWidth; X++)
            {
                int32 Index = X + Y * ClothWidth;

                // Cross springs
                if (X < ClothWidth - 2 && Y < ClothHeight - 2)
                {
                    Springs.Add(FSpring(Index, Index + ClothWidth + 2, FMath::Sqrt(5) * ClothSpacing, 1.3f));
                    Springs.Add(FSpring(Index + 2, Index + ClothWidth, FMath::Sqrt(5) * ClothSpacing, 1.3f));
                }
            }
        }
    }
    else
    {
        // Remove cross springs
        for (int32 s = Springs.Num() - 1; s >= 0; s--)
        {
            FSpring& Spring = Springs[s];
            float ExpectedRestLength1 = FMath::Sqrt(5) * ClothSpacing;
            if (FMath::Abs(Spring.RestLength - ExpectedRestLength1) < 0.1f)
            {
                Springs.RemoveAt(s);
            }
        }
    }
}

void AClothMesh::AdjustClothParameters(int32 NewWidth, int32 NewHeight, int32 NewNumberOfHooks)
{
    ClothWidth = NewWidth;
    ClothHeight = NewHeight;
    NumberOfHooks = NewNumberOfHooks;

    ResetCloth();
}

void AClothMesh::AddCollisionObject(AActor* CollisionActor)
{
    if (CollisionActor)
    {
        CollisionObjects.Add(CollisionActor);
    }
}

void AClothMesh::ApplyGrabForce(FVector GrabPosition, FVector Force)
{
    // Find the closest particle to the grab position
    float MinDistance = FLT_MAX;
    int32 ClosestIndex = -1;

    for (int32 i = 0; i < Particles.Num(); i++)
    {
        float Distance = FVector::Dist(Particles[i].Position, GrabPosition);
        if (Distance < MinDistance)
        {
            MinDistance = Distance;
            ClosestIndex = i;
        }
    }

    if (ClosestIndex != -1)
    {
        // Apply force to the closest particle
        Particles[ClosestIndex].Acceleration += Force;
    }
}

void AClothMesh::OnMouseClick()
{
    APlayerController* PC = GetWorld()->GetFirstPlayerController();
    if (PC)
    {
        FVector WorldLocation, WorldDirection;
        PC->DeprojectMousePositionToWorld(WorldLocation, WorldDirection);

        // Perform a line trace (raycast) to find the closest particle
        FHitResult HitResult;
        FVector End = WorldLocation + (WorldDirection * 10000.0f); // Long enough distance

        // Assuming the cloth mesh has collision enabled
        bool bHit = GetWorld()->LineTraceSingleByChannel(HitResult, WorldLocation, End, ECC_Visibility);
        if (bHit)
        {
            FVector GrabPosition = HitResult.Location;
            FVector Force = WorldDirection * 5000.0f; // Adjust force magnitude as needed

            ApplyGrabForce(GrabPosition, Force);
        }
    }
}

void AClothMesh::EnableTearing(bool bEnable)
{
    bTearingEnabled = bEnable;
}

void AClothMesh::HandleTearing()
{
    // Iterate over springs and remove those that exceed the tear threshold
    for (int32 s = Springs.Num() - 1; s >= 0; s--)
    {
        FSpring& Spring = Springs[s];
        FParticle& pA = Particles[Spring.ParticleA];
        FParticle& pB = Particles[Spring.ParticleB];

        FVector Delta = pB.Position - pA.Position;
        float CurrentLength = Delta.Size();

        if (CurrentLength > TearThreshold * Spring.RestLength)
        {
            Springs.RemoveAt(s);
        }
    }
}

// Curtain Control Functions
void AClothMesh::StartOpening()
{
    if (bIsOpening || bIsClosing)
        return; // Prevent overlapping actions

    if (PinnedParticles.Num() < 2)
    {
        UE_LOG(LogTemp, Warning, TEXT("Not enough hooks to open the curtain."));
        return;
    }

    bIsOpening = true;
    CurrentOpeningHookIndex = PinnedParticles.Num() - 1; // Start with the rightmost hook
    bIsMovingHook = false; // Ready to move the first hook
}

void AClothMesh::StartClosing()
{
    if (bIsOpening || bIsClosing)
        return; // Prevent overlapping actions

    if (PinnedParticles.Num() < 2)
    {
        UE_LOG(LogTemp, Warning, TEXT("Not enough hooks to close the curtain."));
        return;
    }

    bIsClosing = true;
    CurrentClosingHookIndex = PinnedParticles.Num() - 1; // Start with the rightmost hook
    bIsMovingHook = false; // Ready to move the first hook
}

void AClothMesh::SetCurtainTargetOffsets(FVector LeftOffset, FVector RightOffset)
{
    // Ensure no Z-axis offset is applied
    LeftOffset.Z = 0.0f;
    RightOffset.Z = 0.0f;

    CurtainLeftOffset = LeftOffset;
    CurtainRightOffset = RightOffset;

    // Update target positions based on offsets
    TargetHookPositionsLeft.Empty();
    TargetHookPositionsRight.Empty();

    for (int32 i = 0; i < TotalPairs; i++)
    {
        // Left hooks
        FVector OriginalPosLeft = InitialHookPositionsCurtain[i];
        FVector TargetPosLeft = OriginalPosLeft + CurtainLeftOffset;
        // Ensure Z remains unchanged
        TargetPosLeft.Z = OriginalPosLeft.Z;
        TargetHookPositionsLeft.Add(TargetPosLeft);

        // Right hooks
        if (i + TotalPairs < InitialHookPositionsCurtain.Num())
        {
            FVector OriginalPosRight = InitialHookPositionsCurtain[i + TotalPairs];
            FVector TargetPosRight = OriginalPosRight + CurtainRightOffset;
            // Ensure Z remains unchanged
            TargetPosRight.Z = OriginalPosRight.Z;
            TargetHookPositionsRight.Add(TargetPosRight);
        }
    }
}

void AClothMesh::ResetCurtain()
{
    // Reset to initial positions
    for (int32 i = 0; i < PinnedParticles.Num(); i++)
    {
        Particles[PinnedParticles[i]].Position = InitialHookPositionsCurtain[i];
        Particles[PinnedParticles[i]].PreviousPosition = InitialHookPositionsCurtain[i];
    }

    // Update mesh immediately
    UpdateMesh();

    // Clear movement flags
    bIsOpening = false;
    bIsClosing = false;
    bIsMovingHook = false;
    CurrentOpeningHookIndex = 0;
    CurrentClosingHookIndex = 0;
}

void AClothMesh::InitializeCurtainTargets()
{
    FVector DefaultLeftOffset = FVector(-500.0f, 0.0f, 0.0f);
    FVector DefaultRightOffset = FVector(500.0f, 0.0f, 0.0f);

    SetCurtainTargetOffsets(DefaultLeftOffset, DefaultRightOffset);
}

void AClothMesh::UpdateCurtainMovement(float DeltaTime)
{
    if (bIsOpening)
    {
        MoveNextOpeningHook(DeltaTime);
    }

    if (bIsClosing)
    {
        MoveNextClosingHook(DeltaTime);
    }
}

void AClothMesh::MoveNextOpeningHook(float DeltaTime)
{
    if (CurrentOpeningHookIndex < 0)
    {
        bIsOpening = false; // All hooks have been moved
        return;
    }

    if (!bIsMovingHook)
    {
        bIsMovingHook = true; // Start moving the current hook
    }

    // Define target position: Move closer to the previous hook
    FVector CurrentHookInitialPos = InitialHookPositionsCurtain[CurrentOpeningHookIndex];
    FVector PreviousHookPos = (CurrentOpeningHookIndex > 0) ? Particles[PinnedParticles[CurrentOpeningHookIndex - 1]].Position : CurrentHookInitialPos;
    float AdjustableDistance = 50.0f; // Adjustable distance between hooks

    FVector TargetPos = PreviousHookPos + (CurrentHookInitialPos - PreviousHookPos).GetSafeNormal() * AdjustableDistance;

    // Ensure Z remains unchanged
    TargetPos.Z = CurrentHookInitialPos.Z;

    // Move the current hook towards the target position
    Particles[PinnedParticles[CurrentOpeningHookIndex]].Position = FMath::VInterpConstantTo(
        Particles[PinnedParticles[CurrentOpeningHookIndex]].Position,
        TargetPos,
        DeltaTime,
        PairMovementSpeed
    );

    // Check if the hook has reached the target position
    float DistanceToTarget = FVector::Dist(Particles[PinnedParticles[CurrentOpeningHookIndex]].Position, TargetPos);
    if (DistanceToTarget < 1.0f) // Threshold can be adjusted
    {
        // Hook has reached the target
        bIsMovingHook = false;
        CurrentOpeningHookIndex--; // Move to the next hook in the next tick
    }
}

void AClothMesh::MoveNextClosingHook(float DeltaTime)
{
    if (CurrentClosingHookIndex < 0)
    {
        bIsClosing = false; // All hooks have been moved back
        return;
    }

    if (!bIsMovingHook)
    {
        bIsMovingHook = true; // Start moving the current hook
    }

    // Define target position: Move back to the initial position
    FVector TargetPos = InitialHookPositionsCurtain[CurrentClosingHookIndex];

    // Move the current hook towards the target position
    Particles[PinnedParticles[CurrentClosingHookIndex]].Position = FMath::VInterpConstantTo(
        Particles[PinnedParticles[CurrentClosingHookIndex]].Position,
        TargetPos,
        DeltaTime,
        PairMovementSpeed
    );

    // Check if the hook has reached the target position
    float DistanceToTarget = FVector::Dist(Particles[PinnedParticles[CurrentClosingHookIndex]].Position, TargetPos);
    if (DistanceToTarget < 1.0f) // Threshold can be adjusted
    {
        // Hook has reached the target
        bIsMovingHook = false;
        CurrentClosingHookIndex--; // Move to the next hook in the next tick
    }
}


void AClothMesh::InterpolateHookPositions(float DeltaTime, const TArray<FVector>& TargetPositions, bool bIsOpeningFlag)
{
    if (bIsOpeningFlag)
    {
        for (int32 i = 0; i < TargetPositions.Num(); i++)
        {
            if (CurrentOpeningPairIndex >= TotalPairs)
                break;

            int32 LeftHookPairIndex = CurrentOpeningPairIndex;
            int32 RightHookPairIndex = PinnedParticles.Num() - 1 - CurrentOpeningPairIndex;

            // Safety checks
            if (LeftHookPairIndex >= InitialHookPositionsCurtain.Num() || RightHookPairIndex >= InitialHookPositionsCurtain.Num())
            {
                UE_LOG(LogTemp, Error, TEXT("HookPairIndex out of bounds!"));
                continue;
            }

            if (i < TargetPositions.Num())
            {
                FVector InterpolatedTargetPosLeft = TargetPositions[i];
                FVector InterpolatedTargetPosRight = TargetPositions[i];

                // Enforce Z-axis lock
                InterpolatedTargetPosLeft.Z = InitialHookPositionsCurtain[LeftHookPairIndex].Z;
                InterpolatedTargetPosRight.Z = InitialHookPositionsCurtain[RightHookPairIndex].Z;

                // Left hook interpolation
                Particles[PinnedParticles[LeftHookPairIndex]].Position =
                    FMath::VInterpConstantTo(Particles[PinnedParticles[LeftHookPairIndex]].Position, InterpolatedTargetPosLeft, DeltaTime, PairMovementSpeed);

                // Right hook interpolation
                Particles[PinnedParticles[RightHookPairIndex]].Position =
                    FMath::VInterpConstantTo(Particles[PinnedParticles[RightHookPairIndex]].Position, InterpolatedTargetPosRight, DeltaTime, PairMovementSpeed);
            }
        }

        CurrentOpeningPairIndex++;
        if (CurrentOpeningPairIndex >= TargetHookPositionsLeft.Num())
        {
            bIsOpening = false;
        }
    }
    else
    {
        for (int32 i = 0; i < TargetPositions.Num(); i++)
        {
            if (CurrentClosingPairIndex >= TotalPairs)
                break;

            int32 LeftHookPairIndex = CurrentClosingPairIndex;
            int32 RightHookPairIndex = PinnedParticles.Num() - 1 - CurrentClosingPairIndex;

            // Safety checks
            if (LeftHookPairIndex >= InitialHookPositionsCurtain.Num() || RightHookPairIndex >= InitialHookPositionsCurtain.Num())
            {
                UE_LOG(LogTemp, Error, TEXT("HookPairIndex out of bounds!"));
                continue;
            }

            if (i < TargetPositions.Num())
            {
                FVector InterpolatedTargetPosLeft = TargetPositions[i];
                FVector InterpolatedTargetPosRight = TargetPositions[i];

                // Enforce Z-axis lock
                InterpolatedTargetPosLeft.Z = InitialHookPositionsCurtain[LeftHookPairIndex].Z;
                InterpolatedTargetPosRight.Z = InitialHookPositionsCurtain[RightHookPairIndex].Z;

                // Left hook interpolation
                Particles[PinnedParticles[LeftHookPairIndex]].Position =
                    FMath::VInterpConstantTo(Particles[PinnedParticles[LeftHookPairIndex]].Position, InterpolatedTargetPosLeft, DeltaTime, PairMovementSpeed);

                // Right hook interpolation
                Particles[PinnedParticles[RightHookPairIndex]].Position =
                    FMath::VInterpConstantTo(Particles[PinnedParticles[RightHookPairIndex]].Position, InterpolatedTargetPosRight, DeltaTime, PairMovementSpeed);
            }
        }

        CurrentClosingPairIndex++;
        if (CurrentClosingPairIndex >= TargetPositions.Num())
        {
            bIsClosing = false;
        }
    }
}

void AClothMesh::LockTopRowZAxis()
{
    int32 RowsToLock = FMath::Clamp(FMath::CeilToInt(20.0f / ClothSpacing), 1, ClothHeight);

    for (int32 Y = 0; Y < RowsToLock; Y++)
    {
        for (int32 X = 0; X < ClothWidth; X++)
        {
            int32 Index = Y * ClothWidth + X;
            float InitialZ = Particles[Index].Position.Z;
            Particles[Index].Position.Z = InitialZ;
            Particles[Index].PreviousPosition.Z = InitialZ;
        }
    }
}
