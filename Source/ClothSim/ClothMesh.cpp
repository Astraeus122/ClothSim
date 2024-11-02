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
    NumberOfHooks = 6; // Default number of hooks
    GravityStrength = -980.0f; // Corrected to downward direction
    WindStrength = 0.0f;
    WindDirection = FVector(1, 0, 0);
    Damping = 0.99f;
    SelfCollisionRadius = 10.0f;
    SelfCollisionStiffness = 0.5f;
    TearThreshold = 1.5f;
    bTearingEnabled = false;

    // Initialize curtain movement variables
    bIsOpening = false;
    bIsClosing = false;
    CurrentOpeningHookIndex = 0;
    CurrentClosingHookIndex = 0;
    MinHookSpacing = 50.0f;
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
    }
}

void AClothMesh::Tick(float DeltaTime)
{
    Super::Tick(DeltaTime);

    // Simulate the cloth and update the mesh based on particle positions
    SimulateCloth(DeltaTime);
    UpdateMesh();

    // Handle curtain movements based on flags
    if (bIsOpening)
    {
        HandleOpening(DeltaTime);
    }
    if (bIsClosing)
    {
        HandleClosing(DeltaTime);
    }

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

    // Initialize particles above ground
    float InitialZ = 1000.0f; // 10 meters above ground (assuming 1 unit = 1 cm)

    // Ensure a minimum of 3 hooks
    int32 HooksToCreate = FMath::Max(NumberOfHooks, 3);
    int32 MidIndex = ClothWidth / 2;

    // Determine indices for hooks: left end, center, right end, and additional evenly spaced hooks
    TArray<int32> HookIndices;
    HookIndices.Add(0); // Leftmost
    if (HooksToCreate > 2)
    {
        HookIndices.Add(MidIndex); // Center
    }
    HookIndices.Add(ClothWidth - 1); // Rightmost

    int32 HooksBetween = (HooksToCreate - HookIndices.Num()) / 2;
    for (int32 i = 1; i <= HooksBetween; ++i)
    {
        int32 InsertIndexLeft = FMath::Clamp(MidIndex - (i * (MidIndex / (HooksBetween + 1))), 1, ClothWidth - 2);
        HookIndices.Insert(InsertIndexLeft, HookIndices.Num() - i * 2);

        int32 InsertIndexRight = FMath::Clamp(ClothWidth - 1 - (i * (MidIndex / (HooksBetween + 1))), 1, ClothWidth - 2);
        HookIndices.Insert(InsertIndexRight, HookIndices.Num() - i * 2);
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
            }
            else
            {
                Particle.bIsPinned = false;
            }

            Particles.Add(Particle);
        }
    }

    /// Sort PinnedParticles based on their X position (left to right)
    PinnedParticles.Sort([this](const int32 A, const int32 B) -> bool {
        return Particles[A].Position.X < Particles[B].Position.X;
        });

    // Store the sorted initial positions
    InitialHookPositions.Empty();
    for (int32 PinnedIndex : PinnedParticles)
    {
        InitialHookPositions.Add(Particles[PinnedIndex].Position);
    }

    // Initialize curtain target positions
    InitializeCurtainTargets();

    // Create springs (structural, shear, and bend)
    for (int32 Y = 0; Y < ClothHeight; Y++)
    {
        for (int32 X = 0; X < ClothWidth; X++)
        {
            int32 Index = X + Y * ClothWidth;

            // Structural springs
            if (X < ClothWidth - 1)
                Springs.Add(FSpring(Index, Index + 1, ClothSpacing, 2.0f));
            if (Y < ClothHeight - 1)
                Springs.Add(FSpring(Index, Index + ClothWidth, ClothSpacing, 2.0f));

            // Shear springs
            if (X < ClothWidth - 1 && Y < ClothHeight - 1)
            {
                Springs.Add(FSpring(Index, Index + ClothWidth + 1, FMath::Sqrt(2) * ClothSpacing, 1.5f));
                Springs.Add(FSpring(Index + 1, Index + ClothWidth, FMath::Sqrt(2) * ClothSpacing, 1.5f));
            }

            // Bend springs
            if (X < ClothWidth - 2)
                Springs.Add(FSpring(Index, Index + 2, 2 * ClothSpacing, 1.2f));
            if (Y < ClothHeight - 2)
                Springs.Add(FSpring(Index, Index + 2 * ClothWidth, 2 * ClothSpacing, 1.2f));

            // Additional springs for tightly woven structure if enabled
            if (ClothWidth > 3 && ClothHeight > 3)
            {
                // Cross springs
                if (X < ClothWidth - 2 && Y < ClothHeight - 2)
                {
                    Springs.Add(FSpring(Index, Index + ClothWidth + 2, FMath::Sqrt(5) * ClothSpacing, 1.5f));
                    Springs.Add(FSpring(Index + 2, Index + ClothWidth, FMath::Sqrt(5) * ClothSpacing, 1.5f));
                }
            }
        }
    }

    // Initialize movement indices and target positions
    CurrentClosingHookIndex = 0;
    CurrentOpeningHookIndex = 0;
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
            continue;

        FVector Temp = Particles[i].Position;
        // Calculate velocity based on current and previous positions
        Particles[i].Velocity = (Particles[i].Position - Particles[i].PreviousPosition) * Damping;
        // Update position with velocity and acceleration
        Particles[i].Position += Particles[i].Velocity + Particles[i].Acceleration * DeltaTime * DeltaTime;

        // Limit maximum velocity to prevent instability
        float MaxVelocity = 500.0f;
        if (Particles[i].Velocity.Size() > MaxVelocity)
        {
            Particles[i].Velocity = Particles[i].Velocity.GetSafeNormal() * MaxVelocity;
            Particles[i].Position = Particles[i].PreviousPosition + Particles[i].Velocity * DeltaTime;
        }

        Particles[i].PreviousPosition = Temp;
    }

    // Resolve constraints multiple times for stability
    for (int32 Iter = 0; Iter < 15; Iter++) // Increased iterations for better stability
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

            // Clamp the stretch to a maximum of 1.05 times the rest length
            float MaxStretch = 1.05f * Spring.RestLength;
            float ClampedLength = FMath::Min(CurrentLength, MaxStretch);
            float Difference = (ClampedLength - Spring.RestLength) / ClampedLength;

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

    for (FParticle& Particle : Particles)
    {
        // Ground collision
        if (Particle.Position.Z < GroundLevel)
        {
            Particle.Position.Z = GroundLevel;

            // Reflect the Z component of velocity with damping to simulate energy loss
            Particle.Velocity.Z = -Particle.Velocity.Z * 0.5f; // 0.5f is the damping factor
        }

        // Collision with other objects
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
                FVector Correction = Delta.GetSafeNormal() * SelfCollisionRadius;
                Particle.Position = ClosestPoint + Correction;

                // Reflect velocity based on collision normal
                FVector CollisionNormal = Delta.GetSafeNormal();
                Particle.Velocity = Particle.Velocity.MirrorByVector(CollisionNormal) * 0.5f; // Damping factor
            }
        }
    }
}


void AClothMesh::ApplyForces(float DeltaTime)
{
    // Apply gravity
    FVector Gravity = FVector(0, 0, GravityStrength); // Already downward

    // Apply wind with slight oscillation and random variation
    WindTimeAccumulator += DeltaTime;
    float Oscillation = FMath::Sin(WindTimeAccumulator * WindOscillationFrequency);
    FVector DynamicWindDirection = WindDirection.GetSafeNormal();
    // Slightly vary the wind direction over time
    DynamicWindDirection = DynamicWindDirection.RotateAngleAxis(Oscillation * 5.0f, FVector(0, 0, 1));

    // Introduce random variation
    float RandomStrengthVariation = FMath::FRandRange(-10.0f, 10.0f); // Adjust range as needed
    FVector RandomWindDirectionVariation = FVector(FMath::FRandRange(-5.0f, 5.0f), FMath::FRandRange(-5.0f, 5.0f), FMath::FRandRange(-2.0f, 2.0f));
    RandomWindDirectionVariation = RandomWindDirectionVariation.GetClampedToMaxSize(1.0f); // Normalize

    FVector Wind = (DynamicWindDirection + RandomWindDirectionVariation) * (WindStrength + RandomStrengthVariation);

    for (FParticle& Particle : Particles)
    {
        if (Particle.bIsPinned)
            continue;

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
    TMap<int32, TArray<int32>> SpatialHash;

    // Populate the spatial hash
    for (int32 i = 0; i < Particles.Num(); i++)
    {
        int32 Hash = GetCellHash(Particles[i].Position, CellSize);
        SpatialHash.FindOrAdd(Hash).Add(i);
    }

    // Check collisions within each cell and neighboring cells
    for (const auto& Elem : SpatialHash)
    {
        const TArray<int32>& CellParticles = Elem.Value;
        for (int32 i = 0; i < CellParticles.Num(); i++)
        {
            for (int32 j = i + 1; j < CellParticles.Num(); j++)
            {
                int32 IndexA = CellParticles[i];
                int32 IndexB = CellParticles[j];
                FParticle& pA = Particles[IndexA];
                FParticle& pB = Particles[IndexB];

                FVector Delta = pB.Position - pA.Position;
                float Distance = Delta.Size();

                if (Distance < SelfCollisionRadius && Distance > 0.0f)
                {
                    FVector Correction = Delta.GetSafeNormal() * (SelfCollisionRadius - Distance) * 0.8f * SelfCollisionStiffness; // Increased stiffness
                    if (!pA.bIsPinned)
                        pA.Position -= Correction;
                    if (!pB.bIsPinned)
                        pB.Position += Correction;
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
    CurrentClosingHookIndex = PinnedParticles.Num() - 1;
    CurrentOpeningHookIndex = 0;
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

void AClothMesh::HandleClosing(float DeltaTime)
{
    if (CurrentClosingHookIndex < 0)
    {
        bIsClosing = false;
        return;
    }

    for (int32 Offset = 0; Offset < 2; Offset++) // Move two hooks at a time
    {
        if (CurrentClosingHookIndex < 0)
            break;

        int32 PinnedIndex = PinnedParticles[CurrentClosingHookIndex];
        FVector CurrentPos = Particles[PinnedIndex].Position;
        FVector TargetPos = TargetHookPositionsClose[CurrentClosingHookIndex];

        // Allow only X-axis movement for closing, lock Y and Z positions to initial values
        TargetPos.Y = InitialHookPositions[CurrentClosingHookIndex].Y;
        TargetPos.Z = InitialHookPositions[CurrentClosingHookIndex].Z;

        FVector Direction = (TargetPos - CurrentPos).GetSafeNormal();
        float Distance = FVector::Dist(CurrentPos, TargetPos);
        float MoveStep = HookMovementSpeed * DeltaTime;

        if (Distance <= MoveStep)
        {
            Particles[PinnedIndex].Position = TargetPos;
            Particles[PinnedIndex].PreviousPosition = TargetPos;
            CurrentClosingHookIndex--;
        }
        else
        {
            FVector NewPos = CurrentPos + Direction * MoveStep;
            NewPos.Y = InitialHookPositions[CurrentClosingHookIndex].Y;  // Lock Y
            NewPos.Z = InitialHookPositions[CurrentClosingHookIndex].Z;  // Lock Z
            Particles[PinnedIndex].Position = NewPos;
            Particles[PinnedIndex].PreviousPosition = NewPos;
        }
    }
}

void AClothMesh::HandleOpening(float DeltaTime)
{
    if (CurrentOpeningHookIndex >= PinnedParticles.Num())
    {
        bIsOpening = false;
        return;
    }

    for (int32 Offset = 0; Offset < 2; Offset++) // Move two hooks at a time
    {
        if (CurrentOpeningHookIndex >= PinnedParticles.Num())
            break;

        int32 PinnedIndex = PinnedParticles[CurrentOpeningHookIndex];
        FVector CurrentPos = Particles[PinnedIndex].Position;
        FVector TargetPos = TargetHookPositionsOpen[CurrentOpeningHookIndex];

        // Allow only X-axis movement for opening, lock Y and Z positions to initial values
        TargetPos.Y = InitialHookPositions[CurrentOpeningHookIndex].Y;
        TargetPos.Z = InitialHookPositions[CurrentOpeningHookIndex].Z;

        FVector Direction = (TargetPos - CurrentPos).GetSafeNormal();
        float Distance = FVector::Dist(CurrentPos, TargetPos);
        float MoveStep = HookMovementSpeed * DeltaTime;

        if (Distance <= MoveStep)
        {
            Particles[PinnedIndex].Position = TargetPos;
            Particles[PinnedIndex].PreviousPosition = TargetPos;
            CurrentOpeningHookIndex++;
        }
        else
        {
            FVector NewPos = CurrentPos + Direction * MoveStep;
            NewPos.Y = InitialHookPositions[CurrentOpeningHookIndex].Y;  // Lock Y
            NewPos.Z = InitialHookPositions[CurrentOpeningHookIndex].Z;  // Lock Z
            Particles[PinnedIndex].Position = NewPos;
            Particles[PinnedIndex].PreviousPosition = NewPos;
        }
    }
}

void AClothMesh::OpenCurtain()
{
    if (bIsOpening || bIsClosing || PinnedParticles.Num() == 0)
        return;

    bIsOpening = true;
    CurrentOpeningHookIndex = 0; // Start from the leftmost hook
}

void AClothMesh::CloseCurtain()
{
    if (bIsClosing || bIsOpening || PinnedParticles.Num() == 0)
        return;

    bIsClosing = true;
    CurrentClosingHookIndex = PinnedParticles.Num() - 1; // Start from the rightmost hook
}

void AClothMesh::InitializeCurtainTargets()
{
    TargetHookPositionsOpen.Empty();
    TargetHookPositionsClose.Empty();

    float CloseDistance = 200.0f;  // Adjust to control how close hooks come when closed
    float ZOffsetBump = 30.0f;     // Small Z-axis bump for realistic folding
    float YForwardFold = 20.0f;    // Slight depth offset

    // Calculate the leftmost position for all hooks when closing
    FVector LeftmostPosition = InitialHookPositions[0] - FVector(CloseDistance, YForwardFold, ZOffsetBump);

    for (int32 i = 0; i < PinnedParticles.Num(); i++)
    {
        FVector InitialPos = InitialHookPositions[i];

        // Target positions when closing, all moving towards the leftmost side
        FVector TargetPosClose = LeftmostPosition + FVector(0, 0, ZOffsetBump); // Add slight Z bump for realism
        TargetHookPositionsClose.Add(TargetPosClose);

        // Target positions when opening, go back to the initial positions
        TargetHookPositionsOpen.Add(InitialPos);
    }
}
