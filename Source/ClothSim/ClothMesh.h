// ClothMesh.h

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "ProceduralMeshComponent.h"
#include "ClothMesh.generated.h"

// Structure representing a single particle in the cloth simulation
USTRUCT(BlueprintType)
struct FParticle
{
    GENERATED_BODY()

    UPROPERTY(BlueprintReadWrite)
    FVector Position;

    UPROPERTY(BlueprintReadWrite)
    FVector PreviousPosition;

    UPROPERTY(BlueprintReadWrite)
    FVector Acceleration;

    UPROPERTY(BlueprintReadWrite)
    bool bIsPinned;

    UPROPERTY(BlueprintReadWrite)
    FVector Velocity; // Added velocity

    FParticle()
        : Position(FVector::ZeroVector),
        PreviousPosition(FVector::ZeroVector),
        Acceleration(FVector::ZeroVector),
        bIsPinned(false),
        Velocity(FVector::ZeroVector) {}
};

// Structure representing a spring between two particles
USTRUCT(BlueprintType)
struct FSpring
{
    GENERATED_BODY()

    UPROPERTY(BlueprintReadWrite)
    int32 ParticleA;

    UPROPERTY(BlueprintReadWrite)
    int32 ParticleB;

    UPROPERTY(BlueprintReadWrite)
    float RestLength;

    UPROPERTY(BlueprintReadWrite)
    float Stiffness;

    FSpring()
        : ParticleA(0),
        ParticleB(0),
        RestLength(0.f),
        Stiffness(1.f) {}

    FSpring(int32 A, int32 B, float Length, float StiffnessVal)
        : ParticleA(A),
        ParticleB(B),
        RestLength(Length),
        Stiffness(StiffnessVal) {}
};

UCLASS()
class CLOTHSIM_API AClothMesh : public AActor
{
    GENERATED_BODY()

public:
    // Sets default values for this actor's properties
    AClothMesh();

protected:
    // Called when the game starts or when spawned
    virtual void BeginPlay() override;

public:
    // Called every frame
    virtual void Tick(float DeltaTime) override;

    // Cloth Control Functions
    UFUNCTION(BlueprintCallable, Category = "Cloth")
    void ResetCloth();

    UFUNCTION(BlueprintCallable, Category = "Cloth")
    void ReleaseCloth();

    // Additional Features (F1-F6)

    // F1: Fan that creates wind with adjustable speed and direction
    UFUNCTION(BlueprintCallable, Category = "Cloth|Wind")
    void SetWind(FVector NewWindDirection, float NewWindStrength);

    // F2: Toggle tightly woven cloth structure
    UFUNCTION(BlueprintCallable, Category = "Cloth|Structure")
    void ToggleTightlyWoven(bool bEnable);

    // F3: Adjust cloth size and number of hooks
    UFUNCTION(BlueprintCallable, Category = "Cloth|Modularity")
    void AdjustClothParameters(int32 NewWidth, int32 NewHeight, int32 NewNumberOfHooks);

    // F4: Add collision object dynamically
    UFUNCTION(BlueprintCallable, Category = "Cloth|Collision")
    void AddCollisionObject(AActor* CollisionActor);

    // F5: Apply grab force at a specific position
    UFUNCTION(BlueprintCallable, Category = "Cloth|Interaction")
    void ApplyGrabForce(FVector GrabPosition, FVector Force);

    // F6: Enable or disable cloth tearing
    UFUNCTION(BlueprintCallable, Category = "Cloth|Tearing")
    void EnableTearing(bool bEnable);

    // Curtain Control Functions
    /**
     * Starts the process of opening the curtain by pulling it to the side.
     */
    UFUNCTION(BlueprintCallable, Category = "Curtain")
    void StartOpening();

    /**
     * Starts the process of closing the curtain by pulling it back into place.
     */
    UFUNCTION(BlueprintCallable, Category = "Curtain")
    void StartClosing();

    /**
     * Sets the target positions for hooks when the curtain is fully opened.
     * @param LeftOffset The offset to apply to the left hooks.
     * @param RightOffset The offset to apply to the right hooks.
     */
    UFUNCTION(BlueprintCallable, Category = "Curtain")
    void SetCurtainTargetOffsets(FVector LeftOffset, FVector RightOffset);

    /**
     * Resets the curtain to its initial state.
     */
    UFUNCTION(BlueprintCallable, Category = "Curtain")
    void ResetCurtain();

    // Minimum velocity threshold to simulate static friction
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Cloth Physics")
    float StaticFrictionThreshold = 50.0f;

    // Friction coefficients
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Cloth Physics")
    float StaticFrictionCoefficient = 0.5f; 
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Cloth Physics")
    float DynamicFrictionCoefficient = 0.3f;

private:
    // Procedural Mesh Component to render the cloth
    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Cloth", meta = (AllowPrivateAccess = "true"))
    UProceduralMeshComponent* ClothMesh;

    // Function to generate the initial mesh
    void GenerateClothMesh();

    // Function to initialize particles and springs
    void InitializeParticlesAndSprings();

    // Function to simulate cloth physics
    void SimulateCloth(float DeltaTime);

    // Function to apply forces like gravity and wind
    void ApplyForces(float DeltaTime);

    // Store initial positions of hooks to prevent over-stretching
    UPROPERTY()
    TArray<FVector> InitialHookPositions;

    // Function to handle collisions with ground and other objects
    void HandleCollisions();

    // Function to handle self-collisions within the cloth
    void HandleSelfCollisions();

    // Function to update the procedural mesh based on particle positions
    void UpdateMesh();

    // Cloth settings (modifiable via Blueprint or UI)
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Cloth Settings", meta = (AllowPrivateAccess = "true"))
    int32 ClothWidth = 50;

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Cloth Settings", meta = (AllowPrivateAccess = "true"))
    int32 ClothHeight = 40;

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Cloth Settings", meta = (AllowPrivateAccess = "true"))
    float ClothSpacing = 10.0f;

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Cloth Settings", meta = (AllowPrivateAccess = "true"))
    int32 NumberOfHooks = 10;

    // Physics parameters
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Physics Settings", meta = (AllowPrivateAccess = "true"))
    float GravityStrength = -980.0f; // Units: cm/s² (Unreal's default)

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Physics Settings", meta = (AllowPrivateAccess = "true"))
    FVector WindDirection = FVector(1, 0, 0);

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Physics Settings", meta = (AllowPrivateAccess = "true"))
    float WindStrength = 100.0f;

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Physics Settings", meta = (AllowPrivateAccess = "true"))
    float Damping = 0.98f; // Increased damping for stability

    // Self-collision settings
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Collision Settings", meta = (AllowPrivateAccess = "true"))
    float SelfCollisionRadius = 5.0f; // Minimum distance between particles

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Collision Settings", meta = (AllowPrivateAccess = "true"))
    float SelfCollisionStiffness = 0.5f; // Repulsion strength

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Physics Settings", meta = (AllowPrivateAccess = "true"))
    float WindTimeAccumulator = 0.0f;

    // Simulation parameters
    float TimeStep = 0.016f; // Assuming 60 FPS

    // Arrays to hold particles and springs
    UPROPERTY()
    TArray<FParticle> Particles;

    UPROPERTY()
    TArray<FSpring> Springs;

    // Array to hold indices of pinned particles (hooks)
    UPROPERTY()
    TArray<int32> PinnedParticles;

    // Array to hold current vertex positions for the procedural mesh
    UPROPERTY()
    TArray<FVector> VertexPositions;

    // Hook movement parameters
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Hook Settings", meta = (AllowPrivateAccess = "true"))
    float HookMovementSpeed = 100.0f; // Units per second

    // Functions for user interaction (e.g., grabbing the cloth)
    void OnMouseClick();

    // Function to handle tearing
    void HandleTearing();

    // Tearing settings
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Tearing Settings", meta = (AllowPrivateAccess = "true"))
    bool bTearingEnabled = false;

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Tearing Settings", meta = (AllowPrivateAccess = "true"))
    float TearThreshold = 1.5f; // Threshold for tearing springs

    // Collision Objects
    UPROPERTY()
    TArray<AActor*> CollisionObjects;

    // Curtain Control Variables
    // Flags to control hook movement
    bool bMovingLeft;
    bool bMovingRight;

    // Current pair index being moved during closing/opening
    int32 CurrentClosingPairIndex;
    int32 CurrentOpeningPairIndex;

    // Total number of hook pairs
    int32 TotalPairs;

    // Movement speed per pair
    float PairMovementSpeed;

    // Minimum spacing between hooks to prevent overlapping
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Hook Settings", meta = (AllowPrivateAccess = "true"))
    float MinHookSpacing;

    // Target offsets for curtain movement
    FVector CurtainLeftOffset;
    FVector CurtainRightOffset;

    // Initial positions of hooks for curtain movement reset
    TArray<FVector> InitialHookPositionsCurtain;

    // Target positions for hooks when moving left (opening)
    TArray<FVector> TargetHookPositionsLeft;

    // Target positions for hooks when moving right (opening)
    TArray<FVector> TargetHookPositionsRight;

    // Target positions for closing movement (pushing outward)
    TArray<FVector> TargetHookPositionsClosing;

    // Function to update curtain movement each tick
    void UpdateCurtainMovement(float DeltaTime);

    // Function to initialize curtain target positions
    void InitializeCurtainTargets();

    // Function to interpolate hook positions towards target positions
    void InterpolateHookPositions(float DeltaTime, const TArray<FVector>& TargetPositions, bool bIsClosing);

    void ResolveSprings();

    void ApplySpringCorrection(int32 ParticleIndex, const FVector& Correction);

    void ApplyCollisionCorrection(int32 ParticleIndex, const FVector& Correction);

    bool IsTopRowParticle(int32 Index);

    void MoveNextOpeningHook(float DeltaTime);
    void MoveNextClosingHook(float DeltaTime);

    void LockTopRowZAxis();

    // State flags for curtain movement
    bool bIsOpening;
    bool bIsClosing;

    // Current hook index being moved during opening/closing
    int32 CurrentOpeningHookIndex;
    int32 CurrentClosingHookIndex;

    // Flags to indicate if a hook is currently being moved
    bool bIsMovingHook;
};
